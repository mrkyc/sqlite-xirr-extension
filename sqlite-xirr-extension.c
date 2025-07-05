/**
 * @file sqlite-xirr-extension.c
 * @brief SQLite extension for calculating the Internal Rate of Return for a schedule of cash flows
 * that is not necessarily periodic (XIRR).
 *
 * This extension provides XIRR as a user-defined aggregate and window function.
 * It is designed to be efficient for window function usage by using a circular buffer
 * to manage the sliding window of cash flows.
 */
#include <math.h>
#include <sqlite3ext.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

SQLITE_EXTENSION_INIT1

// --- Configuration Constants for XIRR Calculation ---

// The tolerance level for the Newton-Raphson method. This determines the accuracy of the result.
#define XIRR_TOLERANCE 1e-10
// The maximum number of iterations for the Newton-Raphson method to prevent infinite loops.
#define XIRR_MAX_ITERATIONS 100
// The lower bound for the calculated XIRR rate. Prevents rates from becoming excessively negative.
#define XIRR_MIN_RATE -0.9999
// The upper bound for the calculated XIRR rate. Prevents rates from becoming excessively large.
#define XIRR_MAX_RATE 10.0
// An array of default initial guesses (rates) for the Newton-Raphson method to increase the chance of finding a solution.
#define DEFAULT_INITIAL_GUESSES {0.1, 0.05, 0.2, -0.1, 0.01, 0.5}
// The minimum year supported for date parsing, used as the base for date calculations.
#define MIN_DATE_YEAR 1900
// The initial capacity for the dynamic arrays holding cash flow data.
#define INITIAL_CAPACITY 100
// The factor by which the capacity of arrays is increased when they become full.
#define CAPACITY_GROWTH_FACTOR 2

// --- End of Configuration Constants ---

/**
 * @struct XirrData
 * @brief Holds the cash flow data for a window function's state.
 *
 * This structure maintains the state required for calculating XIRR. It uses a
 * circular buffer for efficient addition and removal of cash flows (dates and
 * values), which is crucial for the performance of window functions as the
 * window slides.
 */
typedef struct {
    double *dates;  // Pointer to a dynamic array of cash flow dates (in days).
    double *values; // Pointer to a dynamic array of cash flow values.
    int count;      // The current number of values stored in the buffer.
    int capacity;   // The current allocated capacity of the arrays.
    int head;       // Index of the oldest element (the "front" of the circular buffer).
    int tail;       // Index where the next new element will be inserted (the "back" of the circular buffer).
} XirrData;

/**
 * @struct XirrWindowContext
 * @brief Wrapper structure to hold the context for the XIRR window function.
 *
 * SQLite's aggregate function API uses a context pointer. This structure is
 * what that pointer will point to, holding all the necessary data.
 */
typedef struct {
    XirrData data;                 // The circular buffer and state for cash flows.
    double current_value;          // The last known current value, used as the final cash flow.
    double *custom_starting_rates; // A dynamic array for custom guesses.
    int num_custom_rates;          // Number of custom guesses.
} XirrWindowContext;

// --- Date Helper Functions ---

/**
 * @brief Checks if a given year is a leap year.
 * @param year The year to check.
 * @return 1 if it is a leap year, 0 otherwise.
 */
static int is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

/**
 * @brief Returns the number of days in a given month of a given year.
 * @param year The year.
 * @param month The month (1-12).
 * @return The number of days in the month.
 */
static int days_in_month(int year, int month) {
    int days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month == 2 && is_leap_year(year)) {
        return 29;
    }
    return days[month - 1];
}

/**
 * @brief Converts date components (year, month, day) into a total number of days since MIN_DATE_YEAR.
 * This provides a simple numerical representation for date calculations.
 * @param year The year component of the date.
 * @param month The month component of the date (1-12).
 * @param day The day component of the date.
 * @return The total number of days since the base year (MIN_DATE_YEAR).
 */
static double date_to_days_from_components(int year, int month, int day) {
    long total_days = 0;
    // Accumulate days for full years passed since the base year.
    for (int y = MIN_DATE_YEAR; y < year; y++) {
        total_days += is_leap_year(y) ? 366 : 365;
    }
    // Accumulate days for full months passed in the current year.
    for (int m = 1; m < month; m++) {
        total_days += days_in_month(year, m);
    }
    // Add the remaining days in the current month.
    total_days += day - 1;
    return (double)total_days;
}

/**
 * @brief Parses a date string and converts it to days since MIN_DATE_YEAR.
 * Supports common date formats like "YYYY-MM-DD", "DD/MM/YYYY", and others.
 * @param date_str The date string to parse.
 * @param out_days Pointer to a double where the resulting number of days will be stored.
 * @return SQLITE_OK on successful parsing, SQLITE_MISMATCH on failure.
 */
static int parse_date_string(const char *date_str, double *out_days) {
    if (!date_str)
        return SQLITE_MISMATCH;
    int year, month, day;

    // Try parsing formats like YYYY-MM-DD or YYYY.MM.DD
    if (sscanf(date_str, "%d-%d-%d", &year, &month, &day) == 3 ||
        sscanf(date_str, "%d.%d.%d", &year, &month, &day) == 3) {
        if (year >= MIN_DATE_YEAR && year <= 9999 && month >= 1 && month <= 12 &&
            day >= 1 && day <= days_in_month(year, month)) {
            *out_days = date_to_days_from_components(year, month, day);
            return SQLITE_OK;
        }
    }
    // Try parsing formats like DD/MM/YYYY or MM/DD/YYYY
    if (sscanf(date_str, "%d/%d/%d", &day, &month, &year) == 3 ||
        sscanf(date_str, "%d/%d/%d", &month, &day, &year) == 3) {
        if (year >= MIN_DATE_YEAR && year <= 9999 && month >= 1 && month <= 12 &&
            day >= 1 && day <= days_in_month(year, month)) {
            *out_days = date_to_days_from_components(year, month, day);
            return SQLITE_OK;
        }
    }
    return SQLITE_MISMATCH;
}

// --- Custom Guesses Parser ---

/**
 * @brief Parses a string of custom starting rates separated by a delimiter.
 * @param context The SQLite function context for error reporting.
 * @param rates_str The string containing delimited rates (e.g., "0.1|0.25|-0.1").
 * @param out_rates Pointer to an array of doubles that will be allocated and filled.
 * @param out_count Pointer to an int where the number of parsed rates will be stored.
 * @return SQLITE_OK on success, SQLITE_NOMEM on memory failure, SQLITE_MISMATCH on parsing error.
 */
static int parse_custom_rates(sqlite3_context *context, const char *rates_str, double **out_rates, int *out_count) {
    if (!rates_str || *rates_str == '\0') {
        return SQLITE_MISMATCH;
    }

    // First, count the number of rates to allocate memory.
    int count = 1;
    const char *p = rates_str;
    while (*p) {
        if (*p == '|') {
            count++;
        }
        p++;
    }

    double *rates = (double *)malloc(count * sizeof(double));
    if (!rates) {
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }

    // Now, parse the rates.
    char *str_copy = strdup(rates_str);
    if (!str_copy) {
        free(rates);
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }

    char *token = strtok(str_copy, "|");
    int i = 0;
    while (token != NULL && i < count) {
        char *endptr;
        rates[i] = strtod(token, &endptr);
        if (*endptr != '\0') { // Check if the whole token was a valid double
            free(rates);
            free(str_copy);
            sqlite3_result_error(context, "Invalid number format in initial guesses string", -1);
            return SQLITE_MISMATCH;
        }
        i++;
        token = strtok(NULL, "|");
    }

    free(str_copy);
    *out_rates = rates;
    *out_count = i;
    return SQLITE_OK;
}

// --- XIRR Calculation Helper Functions ---

/**
 * @brief Calculates the Net Present Value (NPV) for a given rate and cash flows.
 * The NPV is the sum of the present values of the individual cash flows.
 * @param rate The discount rate.
 * @param dates Array of cash flow dates (in days).
 * @param values Array of cash flow values.
 * @param count The number of cash flows.
 * @return The calculated NPV, or NAN if inputs are invalid.
 */
static double calculate_npv(double rate, double *dates, double *values, int count) {
    if (count < 2)
        return NAN;
    double npv = 0.0;
    double base_date = dates[0];
    for (int i = 0; i < count; i++) {
        double years = (dates[i] - base_date) / 365.25;
        double factor = pow(1.0 + rate, years);
        if (isinf(factor) || isnan(factor))
            return NAN;
        npv += values[i] / factor;
    }
    return npv;
}

/**
 * @brief Calculates the derivative of the NPV function with respect to the rate.
 * This is required for the Newton-Raphson method.
 * @param rate The discount rate.
 * @param dates Array of cash flow dates (in days).
 * @param values Array of cash flow values.
 * @param count The number of cash flows.
 * @return The calculated derivative of NPV, or NAN if inputs are invalid.
 */
static double calculate_npv_derivative(double rate, double *dates, double *values, int count) {
    if (count < 2)
        return NAN;
    double derivative = 0.0;
    double base_date = dates[0];
    for (int i = 0; i < count; i++) {
        double years = (dates[i] - base_date) / 365.25;
        double factor = pow(1.0 + rate, years + 1.0);
        if (isinf(factor) || isnan(factor))
            return NAN;
        derivative -= years * values[i] / factor;
    }
    return derivative;
}

/**
 * @brief Gets a date/value pair at a logical index in the circular buffer.
 * @param data The XIRR data structure.
 * @param logical_index The 0-based logical index from the start of the window.
 * @param date Pointer to store the retrieved cash flow date.
 * @param value Pointer to store the retrieved cash flow value.
 */
static void get_circular_data(XirrData *data, int logical_index, double *date, double *value) {
    int physical_index = (data->head + logical_index) % data->capacity;
    *date = data->dates[physical_index];
    *value = data->values[physical_index];
}

/**
 * @brief Calculates XIRR using the Newton-Raphson method.
 *
 * This function encapsulates the entire XIRR calculation logic. It creates temporary
 * arrays for cash flows and dates from the context, then iteratively searches for the
 * rate that results in a Net Present Value (NPV) of zero.
 *
 * The function tries multiple starting rates to increase the likelihood of finding a
 * root, as the NPV function can be complex with multiple local extrema. This makes
 * the search for the correct rate more robust.
 *
 * @param ctx The XIRR window context containing all cash flows.
 * @return The calculated XIRR, or NAN if a solution is not found or inputs are invalid.
 */
static double calculate_xirr(XirrWindowContext *ctx) {
    if (!ctx || !ctx->data.values || ctx->data.count < 1) {
        return NAN;
    }

    // Prepare temporary arrays for calculation. The total size is the number of cash flows
    // plus one for the final current_value.
    int total_count = ctx->data.count + 1;
    double *temp_dates = (double *)malloc(total_count * sizeof(double));
    double *temp_values = (double *)malloc(total_count * sizeof(double));

    if (!temp_dates || !temp_values) {
        if (temp_dates)
            free(temp_dates);
        if (temp_values)
            free(temp_values);
        return NAN; // Cannot signal memory error here, so return NAN
    }

    // Copy data from the circular buffer to the temporary linear arrays.
    for (int i = 0; i < ctx->data.count; i++) {
        get_circular_data(&ctx->data, i, &temp_dates[i], &temp_values[i]);
    }

    // Add the final market value as the last cash flow.
    temp_values[ctx->data.count] = ctx->current_value;
    double last_date, last_val;
    get_circular_data(&ctx->data, ctx->data.count - 1, &last_date, &last_val);
    temp_dates[ctx->data.count] = last_date;

    // XIRR requires at least one positive and one negative cash flow.
    int has_positive = 0, has_negative = 0;
    for (int i = 0; i < total_count; i++) {
        if (temp_values[i] > 0)
            has_positive = 1;
        if (temp_values[i] < 0)
            has_negative = 1;
    }

    if (!has_positive || !has_negative) {
        free(temp_values);
        free(temp_dates);
        return NAN;
    }

    // Determine which starting rates to use.
    double default_initial_guesses[] = DEFAULT_INITIAL_GUESSES;
    double *starting_rates = default_initial_guesses;
    int num_starts = sizeof(default_initial_guesses) / sizeof(default_initial_guesses[0]);

    if (ctx->custom_starting_rates && ctx->num_custom_rates > 0) {
        starting_rates = ctx->custom_starting_rates;
        num_starts = ctx->num_custom_rates;
    }

    double result_rate = NAN;

    // Iterate through different starting rates to find a solution.
    for (int start_idx = 0; start_idx < num_starts; start_idx++) {
        double rate = starting_rates[start_idx];
        // Perform Newton-Raphson iteration.
        for (int i = 0; i < XIRR_MAX_ITERATIONS; i++) {
            double npv = calculate_npv(rate, temp_dates, temp_values, total_count);
            double npv_derivative = calculate_npv_derivative(rate, temp_dates, temp_values, total_count);

            if (isnan(npv) || isinf(npv) || isnan(npv_derivative) || isinf(npv_derivative))
                break;

            // If NPV is close to zero, we found the root.
            if (fabs(npv) < XIRR_TOLERANCE) {
                result_rate = rate;
                goto cleanup_and_return;
            }
            // If the derivative is too small, the method is unstable.
            if (fabs(npv_derivative) < XIRR_TOLERANCE)
                break;

            // Newton-Raphson formula.
            double new_rate = rate - npv / npv_derivative;
            if (new_rate < XIRR_MIN_RATE)
                new_rate = XIRR_MIN_RATE;
            if (new_rate > XIRR_MAX_RATE)
                new_rate = XIRR_MAX_RATE;

            // If the rate converges, we found the solution.
            if (fabs(new_rate - rate) < XIRR_TOLERANCE) {
                result_rate = new_rate;
                goto cleanup_and_return;
            }
            rate = new_rate;
        }
    }

cleanup_and_return:
    free(temp_values);
    free(temp_dates);
    return result_rate;
}

// --- Circular Buffer Helper Functions ---

/**
 * @brief Adds a new date/value pair to the end (tail) of the circular buffer.
 * @param data The XIRR data structure.
 * @param date The cash flow date to add.
 * @param value The cash flow value to add.
 */
static void add_to_circular_buffer(XirrData *data, double date, double value) {
    data->dates[data->tail] = date;
    data->values[data->tail] = value;
    data->tail = (data->tail + 1) % data->capacity;
    data->count++;
}

/**
 * @brief Removes a value from the beginning (head) of the circular buffer.
 * @param data The XIRR data structure.
 */
static void remove_from_circular_buffer(XirrData *data) {
    if (data->count == 0)
        return;
    data->head = (data->head + 1) % data->capacity;
    data->count--;
}

// --- Context Management and Result Handling ---

/**
 * @brief Initializes the XirrData structure.
 * @param context The SQLite function context for error reporting.
 * @param data The XirrData structure to initialize.
 * @return SQLITE_OK on success, SQLITE_NOMEM on memory allocation failure.
 */
static int init_xirr_data(sqlite3_context *context, XirrData *data) {
    data->capacity = INITIAL_CAPACITY;
    data->dates = (double *)malloc(data->capacity * sizeof(double));
    data->values = (double *)malloc(data->capacity * sizeof(double));
    data->count = 0;
    data->head = 0;
    data->tail = 0;

    if (!data->dates || !data->values) {
        if (data->dates)
            free(data->dates);
        if (data->values)
            free(data->values);
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }
    return SQLITE_OK;
}

/**
 * @brief Grows the buffers within the XirrData structure.
 * @param context The SQLite function context for error reporting.
 * @param data The XirrData structure to grow.
 * @return SQLITE_OK on success, SQLITE_NOMEM on memory allocation failure.
 */
static int grow_xirr_buffer(sqlite3_context *context, XirrData *data) {
    int new_capacity = data->capacity * CAPACITY_GROWTH_FACTOR;
    double *new_dates = (double *)malloc(new_capacity * sizeof(double));
    double *new_values = (double *)malloc(new_capacity * sizeof(double));

    if (!new_dates || !new_values) {
        if (new_dates)
            free(new_dates);
        if (new_values)
            free(new_values);
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }

    // Copy existing data to the new buffers, maintaining logical order.
    for (int i = 0; i < data->count; i++) {
        get_circular_data(data, i, &new_dates[i], &new_values[i]);
    }

    free(data->dates);
    free(data->values);
    data->dates = new_dates;
    data->values = new_values;
    data->capacity = new_capacity;
    data->head = 0;
    data->tail = data->count;

    return SQLITE_OK;
}

/**
 * @brief Helper to set the result, handling NAN/INF values.
 * @param context The SQLite function context.
 * @param result The double result to set.
 */
static void set_result(sqlite3_context *context, double result) {
    if (isnan(result) || isinf(result)) {
        sqlite3_result_null(context);
    } else {
        sqlite3_result_double(context, result);
    }
}

// --- UNIFIED CALLBACKS FOR AGGREGATE AND WINDOW FUNCTIONS ---

/**
 * @brief The "step" function, called for each row in the aggregate or window frame.
 *
 * This function is called for each row. It adds the new cash flow (value and date)
 * to the XIRR context. It handles context initialization and buffer growth.
 *
 * @param context The SQLite function context.
 * @param argc The number of arguments.
 * @param argv The argument values (date, cashflow, [current_value], [guesses]).
 */
static void xirr_step(sqlite3_context *context, int argc, sqlite3_value **argv) {
    // Validate the number of arguments.
    if (argc < 2 || argc > 4) {
        sqlite3_result_error(context, "XIRR requires 2 to 4 arguments: date, cashflow, [current_value], [guesses]", -1);
        return;
    }

    // Get or create the aggregate context.
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, sizeof(XirrWindowContext));
    if (!ctx) {
        sqlite3_result_error_nomem(context);
        return;
    }

    // Initialize context on the first call.
    if (ctx->data.values == NULL) {
        if (init_xirr_data(context, &ctx->data) != SQLITE_OK)
            return;
        ctx->current_value = 0.0;
        ctx->custom_starting_rates = NULL;
        ctx->num_custom_rates = 0;
    }

    // Ignore rows with NULL required values.
    if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
        return;
    }

    // Grow buffer if it is full.
    if (ctx->data.count >= ctx->data.capacity) {
        if (grow_xirr_buffer(context, &ctx->data) != SQLITE_OK)
            return;
    }

    // --- Argument Parsing ---
    const char *date_str = (const char *)sqlite3_value_text(argv[0]);
    double cashflow = sqlite3_value_double(argv[1]);
    double date_days;

    if (parse_date_string(date_str, &date_days) != SQLITE_OK) {
        sqlite3_result_error(context, "Invalid date format. Use YYYY-MM-DD, DD/MM/YYYY, MM/DD/YYYY, or YYYY.MM.DD", -1);
        return;
    }

    add_to_circular_buffer(&ctx->data, date_days, cashflow);

    // Process optional arguments (3rd and 4th).
    for (int i = 2; i < argc; i++) {
        int arg_type = sqlite3_value_type(argv[i]);
        if (arg_type == SQLITE_NULL)
            continue;

        if (arg_type == SQLITE_FLOAT || arg_type == SQLITE_INTEGER) {
            // This is the current_value
            ctx->current_value = sqlite3_value_double(argv[i]);
        } else if (arg_type == SQLITE_TEXT) {
            // This is the custom guesses string.
            // If custom rates were provided in a previous row for the same aggregate,
            // free the old ones to ensure only the latest set is used.
            if (ctx->custom_starting_rates) {
                free(ctx->custom_starting_rates);
                ctx->custom_starting_rates = NULL;
                ctx->num_custom_rates = 0;
            }
            const char *rates_str = (const char *)sqlite3_value_text(argv[i]);
            if (parse_custom_rates(context, rates_str, &ctx->custom_starting_rates, &ctx->num_custom_rates) != SQLITE_OK) {
                // Error is already set by parse_custom_rates
                return;
            }
        }
    }
}

/**
 * @brief The "inverse" function, called when a row moves out of a window frame.
 * This removes the oldest cash flow from the circular buffer.
 * @param context The SQLite function context.
 * @param argc The number of arguments.
 * @param argv The argument values of the row leaving the window.
 */
static void xirr_inverse(sqlite3_context *context, int argc, sqlite3_value **argv) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (!ctx || !ctx->data.values || ctx->data.count <= 0)
        return;

    remove_from_circular_buffer(&ctx->data);
    if (ctx->data.count == 0) {
        ctx->current_value = 0.0;
    }
}

/**
 * @brief The "value" function, called to get the current result for a window frame.
 * @param context The SQLite function context.
 */
static void xirr_value(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    double result = calculate_xirr(ctx);
    set_result(context, result);
}

/**
 * @brief The "final" function, called to get the final result in aggregate mode.
 *
 * This function calculates the final result and is also responsible for
 * cleaning up any allocated resources for the aggregate context.
 * @param context The SQLite function context.
 */
static void xirr_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    double result = calculate_xirr(ctx);
    set_result(context, result);

    // Clean up allocated memory.
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
    }
}

// --- Extension Initialization ---

/**
 * @brief The main entry point for the SQLite extension.
 *
 * This function is called by SQLite when the extension is loaded. It registers
 * the custom XIRR function (in both lowercase and uppercase forms).
 *
 * @param db The database connection.
 * @param pzErrMsg A pointer to an error message string.
 * @param pApi A pointer to the SQLite API routines.
 * @return SQLITE_OK on success, or an error code on failure.
 */
int sqlite3_extension_init(sqlite3 *db, char **pzErrMsg, const sqlite3_api_routines *pApi) {
    int rc = SQLITE_OK;
    SQLITE_EXTENSION_INIT2(pApi);

    // Register the function for lowercase name "xirr"
    rc = sqlite3_create_window_function(
        db, "xirr", -1, SQLITE_UTF8 | SQLITE_DETERMINISTIC | SQLITE_INNOCUOUS, 0, xirr_step,
        xirr_final, xirr_value, xirr_inverse, NULL);
    if (rc != SQLITE_OK) {
        return rc;
    }

    // Register the function for uppercase name "XIRR"
    rc = sqlite3_create_window_function(
        db, "XIRR", -1, SQLITE_UTF8 | SQLITE_DETERMINISTIC | SQLITE_INNOCUOUS, 0, xirr_step,
        xirr_final, xirr_value, xirr_inverse, NULL);

    return rc;
}
