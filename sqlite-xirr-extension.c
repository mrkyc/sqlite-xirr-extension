/**
 * @file sqlite-xirr-extension.c
 * @brief SQLite extension for calculating the Internal Rate of Return (XIRR) for irregular cash flows.
 *
 * This extension provides four user-defined functions, each available in lowercase and uppercase versions:
 * - xirr(): Calculates a single, stable XIRR value.
 * - xirr_all(): Returns all potential XIRR solutions from a set of initial guesses, separated by pipes.
 * - xirr_unique(): Returns only the unique, valid solutions from xirr_all, separated by pipes.
 * - xirr_scan(): Scans a wide range of rates to find all possible solutions, then returns them sorted and pipe-separated.
 *
 * All functions are implemented as efficient aggregate and window functions. When used as window functions,
 * they can take an optional `start_date` argument. If provided, the functions will return NULL for all rows
 * where the date is earlier than `start_date`, and will only begin performing calculations for rows on or after `start_date`.
 * Importantly, the calculation itself always uses the complete, unfiltered set of cash flows within the current window frame.
 */
#include <ctype.h>
#include <math.h>
#include <sqlite3ext.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

SQLITE_EXTENSION_INIT1

// --- Configuration Constants ---

// The minimum year supported for date parsing, used as the base for date calculations.
#define MIN_DATE_YEAR 1900
// The number of days in a year used for the NPV calculation, accounting for leap years.
#define DAYS_PER_YEAR 365.25
// The tolerance level for the Newton-Raphson method. This determines the accuracy of the result.
#define XIRR_TOLERANCE 1e-10
// A very small number used for derivative tolerance to prevent division by zero, which would halt the Newton-Raphson method.
#define XIRR_DERIVATIVE_TOLERANCE 1e-15
// The maximum number of iterations for the Newton-Raphson method to prevent infinite loops in case of non-convergence.
#define XIRR_MAX_ITERATIONS 100
// The lower bound for the calculated XIRR rate. Prevents rates from becoming excessively negative.
#define XIRR_MIN_RATE -0.9999
// The upper bound for the calculated XIRR rate. Prevents rates from becoming excessively large.
#define XIRR_MAX_RATE 10.0
// An array of default initial guesses (rates) for the Newton-Raphson method to increase the chance of finding a stable solution.
#define DEFAULT_INITIAL_GUESSES {0.1, 0.05, 0.2, -0.1, 0.01, 0.5}
// The default start point used for the XIRR_SCAN calculation.
#define DEFAULT_SCAN_START -0.99
// The default end point used for the XIRR_SCAN calculation.
#define DEFAULT_SCAN_END 10.0
// The default step used for the XIRR_SCAN calculation.
#define DEFAULT_SCAN_STEP 0.01
// The initial capacity for the dynamic arrays holding cash flow data.
#define INITIAL_CAPACITY 100
// The factor by which the capacity of arrays is increased when they become full.
#define CAPACITY_GROWTH_FACTOR 2

// --- Data Structures ---

/**
 * @brief Stores cash flow data in a circular buffer for efficient window function operations.
 *
 * This structure holds cash flow data (dates and values). A circular buffer is used
 * to efficiently add and remove elements from the start and end of the window frame
 * without needing to shift elements in memory.
 * @param dates Pointer to the array of dates (as days since an epoch).
 * @param values Pointer to the array of cash flow values.
 * @param count The current number of items in the buffer.
 * @param capacity The total allocated capacity of the buffer.
 * @param head The index of the first (oldest) item in the circular buffer.
 * @param tail The index where the next new element will be inserted.
 */
typedef struct {
    double *dates, *values;
    size_t count, capacity;
    size_t head, tail;
} XirrData;

/**
 * @brief The context for a single window function invocation, holding all necessary state.
 *
 * This structure holds all the state required for an aggregate or window function call,
 * including the circular buffer for cash flows and any optional user-provided parameters.
 * @param data The circular buffer holding the cash flow data for the current window frame.
 * @param current_value The final cash flow value for the current window frame (often from an optional argument).
 * @param current_row_date The date of the current row being processed by the window function.
 * @param custom_starting_rates An array of user-provided initial guesses for the Newton-Raphson method.
 * @param num_custom_rates The number of custom starting rates.
 * @param scan_config_str A user-provided string to configure the `xirr_scan` function (e.g., "start|end|step").
 * @param calculation_start_date When using a window function, this optional date determines the point from which results are calculated.
 *                               For rows with a date before this, the function returns NULL. The calculation itself, however,
 *                               still uses all cash flows in the current window frame.
 * @param start_date_set A flag (0 or 1) indicating whether a custom `calculation_start_date` has been provided by the user.
 */
typedef struct {
    XirrData data;
    double current_value, current_row_date;
    double *custom_starting_rates;
    size_t num_custom_rates;
    char *scan_config_str;
    double calculation_start_date;
    int start_date_set;
} XirrWindowContext;

/**
 * @brief Represents the result of a core XIRR calculation that can return multiple rates.
 *
 * @param rates A dynamically allocated array of calculated IRR values. If NULL, it indicates an error or no result.
 * @param count The number of rates in the array.
 */
typedef struct {
    double *rates;
    size_t count;
} XirrCalcResult;

/**
 * @struct XirrFunctionGroup
 * @brief Defines a group of related XIRR functions to be registered with the same callbacks.
 */
typedef struct {
    const char **names; // Array of function names/aliases.
    size_t name_count;  // Number of names in the array.
    void (*xStep)(sqlite3_context *, int, sqlite3_value **);
    void (*xFinal)(sqlite3_context *); // Pointer to the xFinal function.
    void (*xValue)(sqlite3_context *); // Pointer to the xValue function.
    void (*xInverse)(sqlite3_context *, int, sqlite3_value **);
} XirrFunctionGroup;

// --- Forward Declarations ---

// Core Calculation Logic
static double *select_starting_rates(const XirrWindowContext *ctx, size_t *out_num_starts);
static double calculate_single_xirr(const XirrWindowContext *ctx);
static XirrCalcResult calculate_all_xirr_solutions(const XirrWindowContext *ctx);
static char *calculate_xirr_scan(sqlite3_context *context, const XirrWindowContext *ctx);

// SQLite Callback Functions
static void xirr_step_unified(sqlite3_context *context, int argc, sqlite3_value **argv);
static void xirr_scan_step(sqlite3_context *context, int argc, sqlite3_value **argv);
static void xirr_inverse_unified(sqlite3_context *context, int argc, sqlite3_value **argv);
static void xirr_final(sqlite3_context *context);
static void xirr_value(sqlite3_context *context);
static void xirr_all_final(sqlite3_context *context);
static void xirr_all_value(sqlite3_context *context);
static void xirr_unique_final(sqlite3_context *context);
static void xirr_unique_value(sqlite3_context *context);
static void xirr_scan_final(sqlite3_context *context);
static void xirr_scan_value(sqlite3_context *context);
static void xirr_destroy(void *pAggregate);

// Helper Function Implementations
static bool has_positive_and_negative(const XirrData *data, double current_value);
static double get_base_date(const XirrData *data, double current_row_date);
static int is_leap_year(int year);
static int days_in_month(int year, int month);
static double date_to_days_from_components(int year, int month, int day);
static bool parse_date_string(const char *date_str, double *out_days);
static bool parse_custom_rates(const char *rates_str, double **out_rates, size_t *out_count);
static bool parse_scan_config(const char *config_str, double *start, double *end, double *step);
static double calculate_npv(double rate, const XirrData *data, double current_row_date, double current_value, double base_date);
static double calculate_npv_derivative(double rate, const XirrData *data, double current_row_date, double current_value, double base_date);
static void get_circular_data(const XirrData *data, size_t logical_index, double *date, double *value);
static void add_to_circular_buffer(XirrData *data, double date, double value);
static void remove_from_circular_buffer(XirrData *data);
static int init_xirr_data(XirrData *data);
static int grow_xirr_buffer(XirrData *data);
static void set_result_double(sqlite3_context *context, double result);
static int compare_doubles_desc(const void *a, const void *b);
static char *format_rates_to_string(double *rates, size_t count, int include_empty, int filter_unique, int sort_desc);
static double find_root_newton_raphson(double start_rate, const XirrData *data, double current_row_date, double current_value, double base_date);

// --- Helper Function Implementations ---

/**
 * @brief Checks if a year is a leap year according to the Gregorian calendar rules.
 * @param year The year to check.
 * @return 1 if it is a leap year, 0 otherwise.
 */
static int is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

/**
 * @brief Returns the number of days in a given month of a specific year.
 * @param year The year, needed to check for leap years in February.
 * @param month The month (1-12).
 * @return The number of days in that month.
 */
static int days_in_month(int year, int month) {
    int d[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return (month == 2 && is_leap_year(year)) ? 29 : d[month - 1];
}

/**
 * @brief Converts a date (year, month, day) to a number of days since the epoch (MIN_DATE_YEAR).
 * @param year Year component.
 * @param month Month component.
 * @param day Day component.
 * @return Total days as a double.
 */
static double date_to_days_from_components(int year, int month, int day) {
    long total_days = 0;
    for (int y = MIN_DATE_YEAR; y < year; y++) {
        total_days += is_leap_year(y) ? 366 : 365;
    }
    for (int m = 1; m < month; m++) {
        total_days += days_in_month(year, m);
    }
    return (double)(total_days + day - 1);
}

/**
 * @brief Parses a date string in "YYYY-MM-DD" format into a count of days since the epoch.
 * @param date_str The date string to parse.
 * @param out_days Pointer to store the resulting number of days.
 * @return `true` on successful parsing, `false` on failure.
 */
static bool parse_date_string(const char *date_str, double *out_days) {
    if (!date_str)
        return false;

    int y, m, d;
    const char *p = date_str;
    char *endptr;

    // Parse year
    y = (int)strtol(p, &endptr, 10);
    if (endptr == p || *endptr != '-')
        return false;
    p = endptr + 1;

    // Parse month
    m = (int)strtol(p, &endptr, 10);
    if (endptr == p || *endptr != '-')
        return false;
    p = endptr + 1;

    // Parse day
    d = (int)strtol(p, &endptr, 10);
    if (endptr == p || *endptr != '\0')
        return false; // Must be end of string

    if (y >= MIN_DATE_YEAR && y <= 9999 && m >= 1 && m <= 12 && d >= 1 && d <= days_in_month(y, m)) {
        *out_days = date_to_days_from_components(y, m, d);
        return true;
    }

    return false;
}

/**
 * @brief Parses a pipe-delimited string of initial guess rates (e.g., "0.1|0.2|-0.05") into an array of doubles.
 * @param rates_str The string to parse.
 * @param out_rates Pointer to the allocated array of rates. The caller is responsible for freeing this.
 * @param out_count Pointer to the number of rates found.
 * @return `true` on success, `false` on failure (e.g., invalid format, memory allocation error).
 */
static bool parse_custom_rates(const char *rates_str, double **out_rates, size_t *out_count) {
    if (!rates_str || *rates_str == '\0') {
        return false;
    }

    size_t count = 0;
    size_t capacity = 8;
    double *rates = (double *)malloc(capacity * sizeof(double));
    if (!rates) {
        *out_rates = NULL;
        return false;
    }

    const char *start = rates_str;
    while (start && *start) {
        char *endptr;
        double rate = strtod(start, &endptr);

        if (endptr == start || (*endptr != '\0' && *endptr != '|')) {
            if (rates) {
                free(rates);
                rates = NULL;
            }
            *out_rates = NULL;
            return false; // Invalid character found
        }

        // Grow the array if needed.
        if (count >= capacity) {
            capacity *= 2;
            double *temp = (double *)realloc(rates, capacity * sizeof(double));
            if (!temp) {
                if (rates) {
                    free(rates);
                    rates = NULL;
                }
                *out_rates = NULL;
                return false; // Memory allocation failed
            }
            rates = temp;
        }

        rates[count++] = rate;

        // Move to the next token.
        start = (*endptr == '|') ? (endptr + 1) : NULL;
    }

    *out_rates = rates;
    *out_count = count;
    return true;
}

/**
 * @brief Parses a pipe-delimited string for scan configuration (`start|end|step`).
 *
 * All parts are optional, allowing for flexible formats where missing values retain their defaults. For example:
 * - `"0.1|5.0|0.05"`: Sets start, end, and step.
 * - `"0.2"`: Sets only start.
 * - `"|10.0"`: Sets only end.
 * - `"||0.01"`: Sets only step.
 *
 * @param config_str The string to parse.
 * @param start Pointer to the start rate (updated on success).
 * @param end Pointer to the end rate (updated on success).
 * @param step Pointer to the step size (updated on success).
 * @return `true` on success, `false` on failure (e.g., trailing invalid characters).
 */
static bool parse_scan_config(const char *config_str, double *start, double *end, double *step) {
    if (!config_str || *config_str == '\0') {
        return true; // Empty string is valid, use defaults.
    }

    const char *p = config_str;
    char *endptr;
    double val;

    // Parse start value.
    val = strtod(p, &endptr);
    if (endptr == p) { // No number found
        // It could be an empty start, e.g., "|10.0"
    } else {
        *start = val;
        p = endptr;
    }

    // Parse end value if present.
    if (*p == '|') {
        p++;
        val = strtod(p, &endptr);
        if (endptr > p) {
            *end = val;
            p = endptr;
        }
    }

    // Parse step value if present.
    if (*p == '|') {
        p++;
        val = strtod(p, &endptr);
        if (endptr > p) {
            *step = val;
            p = endptr;
        }
    }

    // Check for any trailing characters, which would be invalid.
    return *p == '\0';
}

/**
 * @brief Calculates the Net Present Value (NPV) for a given discount rate and series of cash flows.
 * @param rate The discount rate to apply.
 * @param data The circular buffer containing the historical cash flows.
 * @param current_row_date The date of the current row being processed.
 * @param current_value The value of the current row.
 * @param base_date The earliest date in the entire cash flow series, used as the reference for discounting.
 * @return The calculated NPV.
 */
static double calculate_npv(double rate, const XirrData *data, double current_row_date, double current_value, double base_date) {
    double npv = 0.0;

    for (size_t i = 0; i < data->count; i++) {
        double date, value;
        get_circular_data(data, i, &date, &value);
        // Calculate the time difference in years from the first cash flow.
        double years = (date - base_date) / DAYS_PER_YEAR;
        // Calculate the discount factor.
        double factor = pow(1.0 + rate, years);
        if (isinf(factor) || isnan(factor))
            return NAN;
        // Add the discounted cash flow to the total NPV.
        npv += value / factor;
    }

    // Also include the value from the optional `current_value` argument, discounted appropriately.
    double current_years = (current_row_date - base_date) / DAYS_PER_YEAR;
    double current_factor = pow(1.0 + rate, current_years);
    if (isinf(current_factor) || isnan(current_factor))
        return NAN;
    npv += current_value / current_factor;

    return npv;
}

/**
 * @brief Calculates the derivative of the NPV function with respect to the rate.
 *
 * This is required for the Newton-Raphson method, which uses the derivative to find the root of the NPV function.
 * @param rate The discount rate.
 * @param data The circular buffer containing cash flows.
 * @param current_row_date The date of the current row.
 * @param current_value The value of the current row.
 * @param base_date The earliest date in the cash flow series.
 * @return The calculated derivative.
 */
static double calculate_npv_derivative(double rate, const XirrData *data, double current_row_date, double current_value, double base_date) {
    double derivative = 0.0;

    for (size_t i = 0; i < data->count; i++) {
        double date, value;
        get_circular_data(data, i, &date, &value);
        double years = (date - base_date) / DAYS_PER_YEAR;
        double factor = pow(1.0 + rate, years + 1.0);
        if (isinf(factor) || isnan(factor))
            return NAN;
        derivative -= years * value / factor;
    }

    double current_years = (current_row_date - base_date) / DAYS_PER_YEAR;
    double current_factor = pow(1.0 + rate, current_years + 1.0);
    if (isinf(current_factor) || isnan(current_factor))
        return NAN;
    derivative -= current_years * current_value / current_factor;

    return derivative;
}

/**
 * @brief Retrieves a date/value pair from the circular buffer using a logical index.
 * @param data The XirrData struct.
 * @param logical_index The logical index (0 to count-1) to retrieve.
 * @param date Pointer to store the date.
 * @param value Pointer to store the value.
 */
static void get_circular_data(const XirrData *data, size_t logical_index, double *date, double *value) {
    // Map the logical index to the physical index in the circular buffer array.
    size_t phys_idx = (data->head + logical_index) % data->capacity;
    *date = data->dates[phys_idx];
    *value = data->values[phys_idx];
}

/**
 * @brief Adds a new cash flow to the end of the circular buffer.
 * @param data The XirrData struct.
 * @param date The date to add.
 * @param value The value to add.
 */
static void add_to_circular_buffer(XirrData *data, double date, double value) {
    // Add the new data at the tail position.
    data->dates[data->tail] = date;
    data->values[data->tail] = value;
    // Move the tail forward, wrapping around if necessary.
    data->tail = (data->tail + 1) % data->capacity;
    data->count++;
}

/**
 * @brief Removes the oldest cash flow from the front of the circular buffer.
 * @param data The XirrData struct.
 */
static void remove_from_circular_buffer(XirrData *data) {
    if (data->count == 0)
        return;
    // Move the head forward, effectively "removing" the oldest element without shifting memory.
    data->head = (data->head + 1) % data->capacity;
    data->count--;
}

/**
 * @brief Initializes the XirrData structure, allocating initial memory for the buffers.
 * @param data The XirrData struct to initialize.
 * @return SQLITE_OK on success, SQLITE_NOMEM on memory allocation failure.
 */
static int init_xirr_data(XirrData *data) {
    data->capacity = INITIAL_CAPACITY;
    data->count = 0;
    data->head = 0;
    data->tail = 0;
    data->dates = (double *)malloc(data->capacity * sizeof(double));
    data->values = (double *)malloc(data->capacity * sizeof(double));
    if (!data->dates || !data->values) {
        if (data->dates) {
            free(data->dates);
            data->dates = NULL;
        }
        if (data->values) {
            free(data->values);
            data->values = NULL;
        }
        return SQLITE_NOMEM;
    }
    return SQLITE_OK;
}

/**
 * @brief Grows the circular buffer when its capacity is reached.
 *
 * This function allocates a new, larger buffer and copies the existing elements
 * from the old circular layout into a new, linear layout.
 * @param data The XirrData struct to grow.
 * @return SQLITE_OK on success, SQLITE_NOMEM on memory allocation failure.
 */
static int grow_xirr_buffer(XirrData *data) {
    size_t new_capacity = data->capacity * CAPACITY_GROWTH_FACTOR;
    double *new_dates = (double *)malloc(new_capacity * sizeof(double));
    double *new_values = (double *)malloc(new_capacity * sizeof(double));

    if (!new_dates || !new_values) {
        if (new_dates) {
            free(new_dates);
            new_dates = NULL;
        }
        if (new_values) {
            free(new_values);
            new_values = NULL;
        }
        return SQLITE_NOMEM;
    }

    // Copy existing data from the old circular buffer to the new linear arrays.
    // A simple realloc is not sufficient here because the data in the circular
    // buffer might be "wrapped around" (e.g., part at the end of the array
    // and part at the beginning). This loop correctly linearizes the data
    // into the new, larger, contiguous buffer.
    for (size_t i = 0; i < data->count; i++) {
        get_circular_data(data, i, &new_dates[i], &new_values[i]);
    }

    if (data->dates) {
        free(data->dates);
        data->dates = NULL;
    }
    if (data->values) {
        free(data->values);
        data->values = NULL;
    }

    // Update the data pointers and reset head/tail for the new, larger linear layout.
    data->dates = new_dates;
    data->values = new_values;
    data->capacity = new_capacity;
    data->head = 0;
    data->tail = data->count;
    return SQLITE_OK;
}

/**
 * @brief Sets the SQLite result to a double value, handling non-finite cases (NAN, INF) by returning NULL.
 * @param context SQLite context.
 * @param result The double value to set.
 */
static void set_result_double(sqlite3_context *context, double result) {
    if (isnan(result) || isinf(result)) {
        sqlite3_result_null(context);
    } else {
        sqlite3_result_double(context, result);
    }
}

/**
 * @brief Comparison function for qsort to sort doubles in descending order.
 * @param a First value.
 * @param b Second value.
 * @return Comparison result (-1, 0, or 1).
 */
static int compare_doubles_desc(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    if (da < db)
        return 1;
    if (da > db)
        return -1;
    return 0;
}

/**
 * @brief Formats an array of rates into a pipe-delimited string.
 *
 * This function uses a three-pass strategy to safely build the result string:
 * 1. Filter: It iterates through the input rates, filtering them based on the provided
 *    options (uniqueness, validity) and storing the results in a temporary array.
 * 2. Measure: It calculates the exact string length required for the filtered rates.
 * 3. Build: It allocates the precise amount of memory and constructs the final string.
 * This approach avoids buffer overflows and unnecessary, repeated memory allocations.
 *
 * @param rates Array of rates to format.
 * @param count Number of rates in the array.
 * @param include_empty If true, include non-converged rates (NAN) as empty strings between pipes.
 * @param filter_unique If true, filter for unique, valid, non-boundary rates.
 * @param sort_desc If true, sort results in descending order before formatting.
 * @return A dynamically allocated string that must be freed by the caller (or passed to SQLite to manage).
 */
static char *format_rates_to_string(double *rates, size_t count, int include_empty, int filter_unique, int sort_desc) {
    if (!rates || count == 0)
        return NULL;

    if (sort_desc) {
        qsort(rates, count, sizeof(double), compare_doubles_desc);
    }

    // Create a temporary array to hold the rates that will be included in the final string.
    double *final_rates = (double *)malloc(count * sizeof(double));
    if (!final_rates)
        return NULL;
    size_t final_count = 0;

    // First pass: filter the rates according to the function's parameters.
    for (size_t i = 0; i < count; i++) {
        if (isnan(rates[i]) || isinf(rates[i]) || rates[i] <= XIRR_MIN_RATE || rates[i] >= XIRR_MAX_RATE) {
            if (include_empty) {
                final_rates[final_count++] = NAN; // Keep placeholder for empty string
            }
            continue;
        }

        if (filter_unique) {
            int is_duplicate = 0;
            for (size_t j = 0; j < final_count; j++) {
                if (!isnan(final_rates[j]) && fabs(rates[i] - final_rates[j]) < XIRR_TOLERANCE) {
                    is_duplicate = 1;
                    break;
                }
            }
            if (is_duplicate) {
                continue;
            }
        }
        final_rates[final_count++] = rates[i];
    }

    if (final_count == 0) {
        if (final_rates) {
            free(final_rates);
            final_rates = NULL;
        }
        return NULL;
    }

    // Second pass: calculate the exact size needed for the result string.
    size_t total_len = 0;
    for (size_t i = 0; i < final_count; i++) {
        if (!isnan(final_rates[i])) {
            char buffer[64]; // Buffer for a single double-to-string conversion.
            total_len += snprintf(buffer, sizeof(buffer), "%.15g", final_rates[i]);
        }
    }
    // Add space for separators (|).
    total_len += (final_count > 0) ? (final_count - 1) : 0;
    // Add space for the null terminator.
    total_len += 1;

    // Third pass: allocate memory and build the string.
    char *result_string = (char *)malloc(total_len);
    if (!result_string) {
        if (final_rates) {
            free(final_rates);
            final_rates = NULL;
        }
        return NULL;
    }

    char *ptr = result_string;
    size_t remaining_space = total_len;
    for (size_t i = 0; i < final_count; i++) {
        if (i > 0) {
            *ptr++ = '|';
            remaining_space--;
        }
        if (!isnan(final_rates[i])) {
            int written = snprintf(ptr, remaining_space, "%.15g", final_rates[i]);
            if (written > 0) {
                ptr += written;
                remaining_space -= written;
            }
        }
    }
    // Null-terminate the string.
    *ptr = '\0';

    if (final_rates) {
        free(final_rates);
        final_rates = NULL;
    }
    return result_string;
}

/**
 * @brief Runs the Newton-Raphson algorithm to find a precise root of the NPV function (i.e., the XIRR).
 *
 * The method iteratively improves the guess for the rate using the formula:
 *   new_rate = rate - NPV(rate) / NPV'(rate)
 * where NPV' is the derivative of the NPV function. The process stops when the
 * change in the rate is below the tolerance threshold, the NPV is close to zero,
 * or the maximum number of iterations is reached.
 *
 * @param start_rate The initial guess for the rate.
 * @param data The circular buffer containing cash flows.
 * @param current_row_date The date of the current row.
 * @param current_value The value of the current row.
 * @param base_date The earliest date in the series, used as the discounting reference.
 * @return The converged rate, or NAN if the algorithm fails to converge.
 */
static double find_root_newton_raphson(double start_rate, const XirrData *data, double current_row_date, double current_value, double base_date) {
    double rate = start_rate;
    for (int i = 0; i < XIRR_MAX_ITERATIONS; i++) {
        double npv = calculate_npv(rate, data, current_row_date, current_value, base_date);
        double npv_derivative = calculate_npv_derivative(rate, data, current_row_date, current_value, base_date);

        // Stop if the calculation results in non-finite numbers or if the derivative is too close to zero (which would cause division by zero).
        if (isnan(npv) || isinf(npv) || isnan(npv_derivative) || isinf(npv_derivative) || fabs(npv_derivative) < XIRR_DERIVATIVE_TOLERANCE) {
            break;
        }
        // If the NPV is very close to zero, we have found the root.
        if (fabs(npv) < XIRR_TOLERANCE) {
            return rate;
        }

        // Calculate the next guess for the rate using the Newton-Raphson formula.
        double new_rate = rate - npv / npv_derivative;
        // Clamp the new rate to our predefined min/max bounds to prevent divergence.
        if (new_rate < XIRR_MIN_RATE)
            new_rate = XIRR_MIN_RATE;
        if (new_rate > XIRR_MAX_RATE)
            new_rate = XIRR_MAX_RATE;

        // If the change in the rate is negligible, we have converged to a solution.
        if (fabs(new_rate - rate) < XIRR_TOLERANCE) {
            return new_rate;
        }
        rate = new_rate;
    }
    // Return NAN if the algorithm did not converge within the maximum number of iterations.
    return NAN;
}

// --- Core Calculation Logic ---

/**
 * @brief Checks if the cash flow series contains both positive and negative values, a prerequisite for XIRR calculation.
 * @param data The circular buffer of cash flows.
 * @param current_value The final cash flow value for the current window.
 * @return `true` if both positive and negative values are present, `false` otherwise.
 */
static bool has_positive_and_negative(const XirrData *data, double current_value) {
    bool has_pos = current_value > 0;
    bool has_neg = current_value < 0;

    if (has_pos && has_neg)
        return true;

    for (size_t i = 0; i < data->count; i++) {
        double date, value;
        get_circular_data(data, i, &date, &value);
        if (value > 0)
            has_pos = true;
        if (value < 0)
            has_neg = true;
        if (has_pos && has_neg)
            return true;
    }
    return has_pos && has_neg;
}

/**
 * @brief Finds the earliest date in the cash flow series to use as the discounting baseline (t=0).
 * @param data The circular buffer of cash flows.
 * @param current_row_date The date of the current row being processed.
 * @return The earliest date as a double representing days since the epoch.
 */
static double get_base_date(const XirrData *data, double current_row_date) {
    double base_date = current_row_date;
    for (size_t i = 0; i < data->count; i++) {
        double date, value;
        get_circular_data(data, i, &date, &value);
        if (date < base_date) {
            base_date = date;
        }
    }
    return base_date;
}

/**
 * @brief Selects the appropriate array of initial guesses for the Newton-Raphson method.
 *
 * It prioritizes user-provided custom rates. If none are available, it falls back to a
 * predefined set of default guesses to cover various scenarios.
 * @param ctx The window context, which may contain custom rates.
 * @param out_num_starts Pointer to store the number of rates in the returned array.
 * @return A pointer to an array of doubles (initial guesses). This pointer does NOT need to be freed,
 *         as it will point to either a static array or a field within the context.
 */
static double *select_starting_rates(const XirrWindowContext *ctx, size_t *out_num_starts) {
    static double default_initial_guesses[] = DEFAULT_INITIAL_GUESSES;

    if (ctx->custom_starting_rates && ctx->num_custom_rates > 0) {
        *out_num_starts = ctx->num_custom_rates;
        return ctx->custom_starting_rates;
    }

    *out_num_starts = sizeof(default_initial_guesses) / sizeof(default_initial_guesses[0]);
    return default_initial_guesses;
}

/**
 * @brief Finds a single, stable XIRR solution using the Newton-Raphson method.
 *
 * For window functions, if a `calculation_start_date` is set, this function will return NAN
 * for any row with a date before `calculation_start_date`. Otherwise, it calculates the XIRR
 * using all cash flows in the current window.
 * @param ctx The window context with all cash flow data.
 * @return The calculated XIRR, or NAN if a solution is not found or inputs are invalid.
 */
static double calculate_single_xirr(const XirrWindowContext *ctx) {
    // For window functions, return NULL if the current row's date is before the specified start date.
    if (ctx->start_date_set && ctx->current_row_date < ctx->calculation_start_date) {
        return NAN;
    }
    // A valid XIRR requires at least one cash flow and a mix of positive and negative flows.
    if (ctx->data.count == 0 || !has_positive_and_negative(&ctx->data, ctx->current_value)) {
        return NAN;
    }

    size_t num_starts;
    double *starting_rates = select_starting_rates(ctx, &num_starts);
    double base_date = get_base_date(&ctx->data, ctx->current_row_date);

    double result_rate = NAN;
    for (size_t start_idx = 0; start_idx < num_starts; start_idx++) {
        double rate = find_root_newton_raphson(starting_rates[start_idx], &ctx->data, ctx->current_row_date, ctx->current_value, base_date);

        // Accept the first valid, finite rate found.
        if (!isnan(rate) && !isinf(rate) && rate > XIRR_MIN_RATE && rate < XIRR_MAX_RATE) {
            result_rate = rate;
            break; // Found a stable rate, stop searching.
        }
    }

    return result_rate;
}

/**
 * @brief Core engine for `xirr_all` and `xirr_unique`. Calculates XIRR for each starting guess.
 *
 * For window functions, if a `calculation_start_date` is set, this function will return an empty result
 * for any row with a date before `calculation_start_date`. Otherwise, it calculates the XIRR
 * using all cash flows in the current window.
 * @param ctx The window context with all cash flow data.
 * @return A XirrCalcResult struct containing an array of all found rates.
 * @note The caller is responsible for freeing the `result.rates` array.
 */
static XirrCalcResult calculate_all_xirr_solutions(const XirrWindowContext *ctx) {
    XirrCalcResult result = {NULL, 0};

    // For window functions, return an empty result if the current row's date is before the specified start date.
    if (ctx->start_date_set && ctx->current_row_date < ctx->calculation_start_date) {
        return result;
    }

    if (ctx->data.count == 0 || !has_positive_and_negative(&ctx->data, ctx->current_value)) {
        return result;
    }

    size_t num_starts;
    double *starting_rates = select_starting_rates(ctx, &num_starts);
    double base_date = get_base_date(&ctx->data, ctx->current_row_date);

    result.rates = (double *)malloc(num_starts * sizeof(double));
    if (!result.rates) {
        return result; // Allocation failed
    }
    result.count = num_starts;

    for (size_t i = 0; i < num_starts; i++) {
        result.rates[i] = find_root_newton_raphson(starting_rates[i], &ctx->data, ctx->current_row_date, ctx->current_value, base_date);
    }

    return result;
}

/**
 * @brief Core engine for `xirr_scan`. Scans a range of rates to find all sign changes in NPV.
 *
 * For window functions, if a `calculation_start_date` is set, this function will return NULL
 * for any row with a date before `calculation_start_date`. Otherwise, it calculates the XIRR
 * using all cash flows in the current window.
 * @param ctx The window context with all cash flow data.
 * @return A dynamically allocated, pipe-separated string of unique, sorted rates. The caller must free this string.
 */
static char *calculate_xirr_scan(sqlite3_context *context, const XirrWindowContext *ctx) {
    // For window functions, return NULL if the current row's date is before the specified start date.
    if (ctx->start_date_set && ctx->current_row_date < ctx->calculation_start_date) {
        return NULL;
    }

    if (ctx->data.count == 0 || !has_positive_and_negative(&ctx->data, ctx->current_value)) {
        return NULL;
    }

    double scan_start = DEFAULT_SCAN_START, scan_end = DEFAULT_SCAN_END, scan_step = DEFAULT_SCAN_STEP;
    if (ctx->scan_config_str) {
        if (!parse_scan_config(ctx->scan_config_str, &scan_start, &scan_end, &scan_step)) {
            sqlite3_result_error(context, "Invalid characters in scan config string", -1);
            return NULL;
        }
    }

    double *found_rates = NULL;
    size_t found_count = 0, found_capacity = 0;

    double base_date = get_base_date(&ctx->data, ctx->current_row_date);
    double prev_npv = calculate_npv(scan_start, &ctx->data, ctx->current_row_date, ctx->current_value, base_date);

    // Scan a range of interest rates to find potential IRR solutions.
    // The strategy is to detect a sign change in the Net Present Value (NPV)
    // between two consecutive steps. A sign change (e.g., from + to -) implies
    // that a root of the NPV function, which is an IRR, exists between those two points.
    for (double r = scan_start + scan_step; r <= scan_end; r += scan_step) {
        double current_npv = calculate_npv(r, &ctx->data, ctx->current_row_date, ctx->current_value, base_date);
        if (prev_npv * current_npv < 0) {
            // Use Newton-Raphson to precisely locate the root within this small interval.
            // We start the search from the midpoint of the interval where the sign change occurred.
            double root = find_root_newton_raphson(r - (scan_step / 2.0), &ctx->data, ctx->current_row_date, ctx->current_value, base_date);

            // Check if the found root is a valid, finite number within our defined bounds.
            if (!isnan(root) && !isinf(root) && root > XIRR_MIN_RATE && root < XIRR_MAX_RATE) {
                // Before adding the root, ensure it's not a duplicate of one we've already found.
                int is_duplicate = 0;
                for (size_t i = 0; i < found_count; i++) {
                    if (fabs(root - found_rates[i]) < XIRR_TOLERANCE) {
                        is_duplicate = 1;
                        break;
                    }
                }
                if (!is_duplicate) {
                    // If the array is full, grow it.
                    if (found_count >= found_capacity) {
                        found_capacity = (found_capacity == 0) ? 10 : found_capacity * 2;
                        double *temp = realloc(found_rates, found_capacity * sizeof(double));
                        if (!temp) {
                            if (found_rates) {
                                free(found_rates);
                                found_rates = NULL;
                            }
                            sqlite3_result_error_nomem(context);
                            return NULL;
                        }
                        found_rates = temp;
                    }
                    // Add the new, unique root to our list.
                    found_rates[found_count++] = root;
                }
            }
        }
        // The current NPV becomes the previous NPV for the next iteration.
        prev_npv = current_npv;
    }

    char *result_string = format_rates_to_string(found_rates, found_count, 0, 0, 1);
    if (found_rates) {
        free(found_rates);
        found_rates = NULL;
    }
    return result_string;
}

// --- SQLite Callback Functions ---
/**
 * @brief Initializes the XirrWindowContext structure upon first use.
 * @param ctx Pointer to the XirrWindowContext structure to initialize.
 */
static void init_xirr_window_context(XirrWindowContext *ctx) {
    ctx->data.values = NULL;
    ctx->data.dates = NULL;
    ctx->data.count = 0;
    ctx->data.capacity = 0;
    ctx->data.head = 0;
    ctx->data.tail = 0;
    ctx->current_value = 0.0;
    ctx->current_row_date = 0.0;
    ctx->custom_starting_rates = NULL;
    ctx->num_custom_rates = 0;
    ctx->scan_config_str = NULL;
    ctx->calculation_start_date = 0.0;
    ctx->start_date_set = 0;
}

/**
 * @brief Common processing logic for the xStep function of all XIRR variants.
 *
 * This function handles context initialization, memory allocation, and parsing of the
 * core date and value arguments, adding them to the circular buffer.
 * @param context The SQLite function context.
 * @param argc The number of arguments.
 * @param argv The array of argument values.
 * @return SQLITE_OK on success, SQLITE_NOMEM on memory error, SQLITE_ERROR on other failures.
 */
static int xirr_step_common(sqlite3_context *context, int argc, sqlite3_value **argv) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, sizeof(XirrWindowContext));
    if (!ctx) {
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }

    // Initialize the context on the first call for this window.
    if (ctx->data.values == NULL) {
        init_xirr_window_context(ctx);
        if (init_xirr_data(&ctx->data) != SQLITE_OK) {
            sqlite3_result_error_nomem(context);
            return SQLITE_NOMEM;
        }
    }

    if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
        return SQLITE_OK; // Ignore NULL inputs gracefully.
    }

    // --- Argument Validation ---
    // Validate arg 1: Date
    if (sqlite3_value_type(argv[0]) != SQLITE_TEXT) {
        sqlite3_result_error(context, "Argument 1 (date) must be a text string in 'YYYY-MM-DD' format.", -1);
        return SQLITE_ERROR;
    }
    double date_days;
    const char *date_str = (const char *)sqlite3_value_text(argv[0]);
    if (!parse_date_string(date_str, &date_days)) {
        char err_msg[128];
        snprintf(err_msg, sizeof(err_msg), "Argument 1 (date) has an invalid format. Expected 'YYYY-MM-DD', but got '%s'.", date_str);
        sqlite3_result_error(context, err_msg, -1);
        return SQLITE_ERROR;
    }

    // Validate arg 2: Value
    int value_type = sqlite3_value_type(argv[1]);
    if (value_type != SQLITE_FLOAT && value_type != SQLITE_INTEGER) {
        sqlite3_result_error(context, "Argument 2 (value) must be a number.", -1);
        return SQLITE_ERROR;
    }
    double value = sqlite3_value_double(argv[1]);
    // --- End Argument Validation ---

    // Grow the circular buffer if it's full.
    if (ctx->data.count >= ctx->data.capacity) {
        if (grow_xirr_buffer(&ctx->data) != SQLITE_OK) {
            sqlite3_result_error_nomem(context);
            return SQLITE_NOMEM;
        }
    }

    // Add the cash flow to the buffer and set current row date.
    add_to_circular_buffer(&ctx->data, date_days, value);
    ctx->current_row_date = date_days; // Set for window function logic

    return SQLITE_OK;
}

/**
 * @brief Unified step function for `xirr`, `xirr_all`, `xirr_unique`.
 *
 * This function is called for each row in the aggregate/window. It parses all arguments,
 * including optional ones, and updates the window context.
 * @param context The SQLite function context.
 * @param argc The number of arguments passed to the function.
 * @param argv The array of argument values.
 */
static void xirr_step_unified(sqlite3_context *context, int argc, sqlite3_value **argv) {
    if (argc < 2 || argc > 5) {
        sqlite3_result_error(context, "XIRR functions require 2 to 5 arguments: xirr(date, value, [current_value], [start_date], [custom_guesses]).", -1);
        return;
    }

    if (xirr_step_common(context, argc, argv) != SQLITE_OK) {
        return;
    }

    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx)
        return;

    // --- Optional Argument Validation ---
    bool has_numeric_opt = false;
    bool has_start_date_opt = false;
    bool has_guesses_opt = false;

    for (int i = 2; i < argc; i++) {
        int arg_type = sqlite3_value_type(argv[i]);
        if (arg_type == SQLITE_NULL)
            continue;

        if (arg_type == SQLITE_FLOAT || arg_type == SQLITE_INTEGER) {
            if (has_numeric_opt) {
                sqlite3_result_error(context, "Provide at most one numeric optional argument (for current_value).", -1);
                return;
            }
            ctx->current_value = sqlite3_value_double(argv[i]);
            has_numeric_opt = true;
        } else if (arg_type == SQLITE_TEXT) {
            const char *text = (const char *)sqlite3_value_text(argv[i]);
            double date_days_opt;

            if (parse_date_string(text, &date_days_opt)) {
                if (has_start_date_opt) {
                    sqlite3_result_error(context, "Provide at most one date string optional argument (for start_date).", -1);
                    return;
                }
                ctx->calculation_start_date = date_days_opt;
                ctx->start_date_set = 1;
                has_start_date_opt = true;
            } else {
                if (has_guesses_opt) {
                    sqlite3_result_error(context, "Provide at most one non-date text optional argument (for custom_guesses).", -1);
                    return;
                }
                double *new_rates = NULL;
                size_t new_num_rates = 0;
                if (!parse_custom_rates(text, &new_rates, &new_num_rates)) {
                    sqlite3_result_error(context, "Invalid format for custom guesses string. Use pipe-separated numbers (e.g., '0.1|0.2').", -1);
                    return;
                }
                if (ctx->custom_starting_rates) {
                    free(ctx->custom_starting_rates);
                    ctx->custom_starting_rates = NULL;
                }
                ctx->custom_starting_rates = new_rates;
                ctx->num_custom_rates = new_num_rates;
                has_guesses_opt = true;
            }
        }
    }
}

/**
 * @brief Step function for `xirr_scan`.
 *
 * Similar to `xirr_step_unified`, but it handles the `scan_config` string instead of custom guesses.
 * @param context The SQLite function context.
 * @param argc The number of arguments passed to the function.
 * @param argv The array of argument values.
 */
static void xirr_scan_step(sqlite3_context *context, int argc, sqlite3_value **argv) {
    if (argc < 2 || argc > 5) {
        sqlite3_result_error(context, "xirr_scan requires 2 to 5 arguments: xirr_scan(date, value, [current_value], [start_date], [scan_config]).", -1);
        return;
    }

    if (xirr_step_common(context, argc, argv) != SQLITE_OK) {
        return;
    }

    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx)
        return;

    // --- Optional Argument Validation ---
    bool has_numeric_opt = false;
    bool has_start_date_opt = false;
    bool has_scan_config_opt = false;

    for (int i = 2; i < argc; i++) {
        int arg_type = sqlite3_value_type(argv[i]);
        if (arg_type == SQLITE_NULL)
            continue;

        if (arg_type == SQLITE_FLOAT || arg_type == SQLITE_INTEGER) {
            if (has_numeric_opt) {
                sqlite3_result_error(context, "Provide at most one numeric optional argument (for current_value).", -1);
                return;
            }
            ctx->current_value = sqlite3_value_double(argv[i]);
            has_numeric_opt = true;
        } else if (arg_type == SQLITE_TEXT) {
            const char *text = (const char *)sqlite3_value_text(argv[i]);
            double date_days;

            if (parse_date_string(text, &date_days)) {
                if (has_start_date_opt) {
                    sqlite3_result_error(context, "Provide at most one date string optional argument (for start_date).", -1);
                    return;
                }
                ctx->calculation_start_date = date_days;
                ctx->start_date_set = 1;
                has_start_date_opt = true;
            } else {
                if (has_scan_config_opt) {
                    sqlite3_result_error(context, "Provide at most one non-date text optional argument (for scan_config).", -1);
                    return;
                }
                if (ctx->scan_config_str) {
                    free(ctx->scan_config_str);
                    ctx->scan_config_str = NULL;
                }
                size_t len = strlen(text) + 1;
                ctx->scan_config_str = (char *)malloc(len);
                if (ctx->scan_config_str) {
                    memcpy(ctx->scan_config_str, text, len);
                } else {
                    sqlite3_result_error_nomem(context);
                    return;
                }
                has_scan_config_opt = true;
            }
        }
    }
}

/**
 * @brief Unified inverse function for all XIRR variants (for window functions).
 *
 * This is called when a row moves out of the window frame. It removes the oldest
 * cash flow from the circular buffer to keep the window's data current.
 * @param context The SQLite function context.
 * @param argc The number of arguments passed to the function.
 * @param argv The array of argument values for the row leaving the window.
 */
static void xirr_inverse_unified(sqlite3_context *context, int argc, sqlite3_value **argv) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (ctx && ctx->data.values && ctx->data.count > 0) {
        remove_from_circular_buffer(&ctx->data);
        if (ctx->data.count == 0) {
            ctx->current_value = 0.0;
        }
    }
}

/**
 * @brief Final callback for `xirr` (aggregate mode).
 *
 * This function is called once at the end of a group in an aggregate query (e.g., with GROUP BY)
 * to compute and return the final result for that group.
 * @param context The SQLite function context.
 */
static void xirr_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx) {
        return;
    }
    double result = calculate_single_xirr(ctx);
    set_result_double(context, result);
}

/**
 * @brief Value callback for `xirr` (window function mode).
 *
 * This function is called for each row within a window partition to calculate and return
 * the result for the current state of the window frame.
 * @param context The SQLite function context.
 */
static void xirr_value(sqlite3_context *context) {
    xirr_final(context);
}

/**
 * @brief Final callback for `xirr_all` (aggregate mode).
 *
 * This function is called once at the end of a group in an aggregate query (e.g., with GROUP BY)
 * to compute and return the final pipe-separated string of all found rates.
 * @param context The SQLite function context.
 */
static void xirr_all_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx) {
        return;
    }
    XirrCalcResult result = calculate_all_xirr_solutions(ctx);

    if (!result.rates) {
        sqlite3_result_null(context);
        return;
    }

    char *result_str = format_rates_to_string(result.rates, result.count, 0, 1, 0);

    if (result_str) {
        sqlite3_result_text(context, result_str, -1, free);
    } else {
        sqlite3_result_null(context);
    }

    if (result.rates) {
        free(result.rates);
        result.rates = NULL;
    }
}

/**
 * @brief Value callback for `xirr_all` (window function mode).
 *
 * This function is called for each row within a window partition to calculate and return
 * the result for the current state of the window frame. It reuses the finalization logic.
 * @param context The SQLite function context.
 */
static void xirr_all_value(sqlite3_context *context) {
    xirr_all_final(context);
}

/**
 * @brief Final callback for `xirr_unique` (aggregate mode).
 *
 * This function is called once at the end of a group in an aggregate query (e.g., with GROUP BY)
 * to compute and return the final pipe-separated string of unique, valid rates.
 * @param context The SQLite function context.
 */
static void xirr_unique_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx) {
        return;
    }
    XirrCalcResult result = calculate_all_xirr_solutions(ctx);

    if (!result.rates) {
        sqlite3_result_null(context);
        return;
    }

    char *result_str = format_rates_to_string(result.rates, result.count, 0, 1, 0);

    if (result_str) {
        sqlite3_result_text(context, result_str, -1, free);
    } else {
        sqlite3_result_null(context);
    }

    if (result.rates) {
        free(result.rates);
        result.rates = NULL;
    }
}

/**
 * @brief Value callback for `xirr_unique` (window function mode).
 *
 * This function is called for each row within a window partition to calculate and return
 * the result for the current state of the window frame. It reuses the finalization logic.
 * @param context The SQLite function context.
 */
static void xirr_unique_value(sqlite3_context *context) {
    xirr_unique_final(context);
}

/**
 * @brief Final callback for `xirr_scan` (aggregate mode).
 *
 * This function is called once at the end of a group in an aggregate query (e.g., with GROUP BY)
 * to compute and return the final pipe-separated string of all rates found by scanning.
 * @param context The SQLite function context.
 */
static void xirr_scan_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx) {
        return;
    }
    char *result_str = calculate_xirr_scan(context, ctx);

    if (result_str) {
        sqlite3_result_text(context, result_str, -1, free);
    } else {
        sqlite3_result_null(context);
    }
}

/**
 * @brief Value callback for `xirr_scan` (window function mode).
 *
 * This function is called for each row within a window partition to calculate and return
 * the result for the current state of the window frame. It reuses the finalization logic.
 * @param context The SQLite function context.
 */
static void xirr_scan_value(sqlite3_context *context) {
    xirr_scan_final(context);
}

/**
 * @brief Destructor for the aggregate context.
 *
 * This function is registered with SQLite and is guaranteed to be called
 * when the aggregation is complete, even if the query is aborted or encounters an error.
 * It ensures that all dynamically allocated memory within the context is freed, preventing memory leaks.
 * @param pAggregate The aggregate context to be destroyed.
 */
static void xirr_destroy(void *pAggregate) {
    XirrWindowContext *ctx = (XirrWindowContext *)pAggregate;
    if (ctx) {
        if (ctx->data.dates) {
            free(ctx->data.dates);
            ctx->data.dates = NULL;
        }
        if (ctx->data.values) {
            free(ctx->data.values);
            ctx->data.values = NULL;
        }
        if (ctx->custom_starting_rates) {
            free(ctx->custom_starting_rates);
            ctx->custom_starting_rates = NULL;
        }
        if (ctx->scan_config_str) {
            free(ctx->scan_config_str);
            ctx->scan_config_str = NULL;
        }
    }
}

/**
 * @brief Helper function to register a unified XIRR function group (lowercase and uppercase).
 * @param db The database connection.
 * @param group The function group to register.
 * @return SQLITE_OK on success, or an error code on failure.
 */
static int register_xirr_function_group(sqlite3 *db, const XirrFunctionGroup *group) {
    int rc = SQLITE_OK;
    int flags = SQLITE_UTF8 | SQLITE_DETERMINISTIC | SQLITE_INNOCUOUS;

    for (size_t i = 0; i < group->name_count; i++) {
        const char *name = group->names[i];
        rc = sqlite3_create_window_function(db, name, -1, flags, 0, group->xStep, group->xFinal, group->xValue, group->xInverse, xirr_destroy);
        if (rc != SQLITE_OK)
            return rc;

        // Create and register the uppercase version for case-insensitivity.
        size_t name_len = strlen(name);
        char *upper_name = malloc(name_len + 1);
        if (!upper_name)
            return SQLITE_NOMEM;
        for (size_t j = 0; j < name_len; j++) {
            upper_name[j] = toupper((unsigned char)name[j]);
        }
        upper_name[name_len] = '\0';

        rc = sqlite3_create_window_function(db, upper_name, -1, flags, 0, group->xStep, group->xFinal, group->xValue, group->xInverse, xirr_destroy);
        if (upper_name) {
            free(upper_name);
            upper_name = NULL;
        }
        if (rc != SQLITE_OK)
            return rc;
    }
    return SQLITE_OK;
}

/**
 * @brief Main entry point for the SQLite extension. Registers all XIRR functions.
 * @param db The database connection handle.
 * @param pzErrMsg A pointer to an error message string, to be set by SQLite on error.
 * @param pApi A pointer to the SQLite API routines.
 * @return SQLITE_OK on success, or an error code on failure.
 */
int sqlite3_extension_init(sqlite3 *db, char **pzErrMsg, const sqlite3_api_routines *pApi) {
    int rc = SQLITE_OK;
    SQLITE_EXTENSION_INIT2(pApi);

    const char *xirr_names[] = {"xirr"};
    const char *xirr_all_names[] = {"xirr_all"};
    const char *xirr_unique_names[] = {"xirr_unique"};
    const char *xirr_scan_names[] = {"xirr_scan"};

    XirrFunctionGroup functions_to_register[] = {
        {xirr_names, sizeof(xirr_names) / sizeof(xirr_names[0]), xirr_step_unified, xirr_final, xirr_value, xirr_inverse_unified},
        {xirr_all_names, sizeof(xirr_all_names) / sizeof(xirr_all_names[0]), xirr_step_unified, xirr_all_final, xirr_all_value, xirr_inverse_unified},
        {xirr_unique_names, sizeof(xirr_unique_names) / sizeof(xirr_unique_names[0]), xirr_step_unified, xirr_unique_final, xirr_unique_value, xirr_inverse_unified},
        {xirr_scan_names, sizeof(xirr_scan_names) / sizeof(xirr_scan_names[0]), xirr_scan_step, xirr_scan_final, xirr_scan_value, xirr_inverse_unified}};

    size_t num_groups = sizeof(functions_to_register) / sizeof(functions_to_register[0]);
    for (size_t i = 0; i < num_groups; i++) {
        rc = register_xirr_function_group(db, &functions_to_register[i]);
        if (rc != SQLITE_OK) {
            return rc;
        }
    }

    return rc;
}