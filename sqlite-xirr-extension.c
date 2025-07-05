/**
 * @file sqlite-xirr-extension.c
 * @brief SQLite extension for calculating the Internal Rate of Return (XIRR) for irregular cash flows.
 *
 * This extension provides four user-defined functions:
 * - xirr(): Calculates a single, stable XIRR value.
 * - xirr_all(): Returns all potential XIRR solutions from a set of initial guesses.
 * - xirr_unique(): Returns only the unique, non-boundary solutions from xirr_all.
 * - xirr_scan(): Scans a wide range of rates to find all possible solutions, then returns them sorted.
 *
 * All functions are implemented as efficient aggregate and window functions.
 */
#include <math.h>
#include <sqlite3ext.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

SQLITE_EXTENSION_INIT1

// --- Configuration Constants ---

// The minimum year supported for date parsing, used as the base for date calculations.
#define MIN_DATE_YEAR 1900
// The number of days in a year used for the NPV calculation.
#define DAYS_PER_YEAR 365.25
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
 * @brief Stores cash flow data in a circular buffer.
 *
 * A circular buffer is an efficient data structure for window functions,
 * where items are added to one end (tail) and removed from the other (head)
 * as the window slides over the data.
 * @param dates Pointer to the array of dates (as days since an epoch).
 * @param values Pointer to the array of cash flow values.
 * @param count The current number of items in the buffer.
 * @param capacity The total capacity of the buffer.
 * @param head The index of the first (oldest) element.
 * @param tail The index where the next element will be inserted.
 */
typedef struct {
    double *dates, *values;
    int count, capacity, head, tail;
} XirrData;

/**
 * @brief The context for a single window function invocation.
 *
 * This structure holds all the state required for an aggregate or window function call,
 * including the circular buffer for cash flows and any optional parameters.
 * @param data The circular buffer holding the cash flow data.
 * @param current_value The final cash flow value for the current window frame.
 * @param custom_starting_rates An array of user-provided initial guesses for the Newton-Raphson method.
 * @param num_custom_rates The number of custom starting rates.
 * @param scan_config_str A user-provided string to configure the `xirr_scan` function (e.g., "start|end|step").
 */
typedef struct {
    XirrData data;
    double current_value, *custom_starting_rates;
    int num_custom_rates;
    char *scan_config_str;
} XirrWindowContext;

/**
 * @brief Represents the result of a core XIRR calculation.
 *
 * @param rates A dynamically allocated array of calculated IRR values.
 * @param count The number of rates in the array.
 */
typedef struct {
    double *rates;
    int count;
} XirrCalcResult;

// --- Forward Declarations ---

static void xirr_destroy(void *pAggregate);
static int prepare_temp_cashflow_data(XirrWindowContext *ctx, double **out_dates, double **out_values, int *out_count);
static XirrCalcResult calculate_xirr_core(XirrWindowContext *ctx);
static char *calculate_xirr_scan(XirrWindowContext *ctx);
static double find_root_newton_raphson(double start_rate, double *dates, double *values, int count);
static void xirr_step_unified(sqlite3_context *context, int argc, sqlite3_value **argv);
static void xirr_inverse_unified(sqlite3_context *context, int argc, sqlite3_value **argv);
static void xirr_final(sqlite3_context *context);
static void xirr_value(sqlite3_context *context);
static void xirr_all_final(sqlite3_context *context);
static void xirr_all_value(sqlite3_context *context);
static void xirr_unique_final(sqlite3_context *context);
static void xirr_unique_value(sqlite3_context *context);
static void xirr_scan_step(sqlite3_context *context, int argc, sqlite3_value **argv);
static void xirr_scan_final(sqlite3_context *context);
static void xirr_scan_value(sqlite3_context *context);
static int is_leap_year(int year);
static int days_in_month(int year, int month);
static double date_to_days_from_components(int year, int month, int day);
static int parse_date_string(const char *date_str, double *out_days);
static int parse_custom_rates(sqlite3_context *context, const char *rates_str, double **out_rates, int *out_count);
static int parse_scan_config(sqlite3_context *context, const char *config_str, double *start, double *end, double *step);
static double calculate_npv(double rate, double *dates, double *values, int count);
static double calculate_npv_derivative(double rate, double *dates, double *values, int count);
static void get_circular_data(XirrData *data, int logical_index, double *date, double *value);
static void add_to_circular_buffer(XirrData *data, double date, double value);
static void remove_from_circular_buffer(XirrData *data);
static int init_xirr_data(sqlite3_context *context, XirrData *data);
static int grow_xirr_buffer(sqlite3_context *context, XirrData *data);
static void set_result_double(sqlite3_context *context, double result);
static int compare_doubles_desc(const void *a, const void *b);
static char *format_rates_to_string(double *rates, int count, int include_empty, int filter_unique, int sort_desc);

// --- Core Calculation Logic ---

/**
 * @brief Core engine for `xirr`, `xirr_all`, `xirr_unique`.
 * @param ctx The window context with all cash flow data.
 * @return A XirrCalcResult struct containing an array of all found rates (or NAN for non-convergence).
 * @note The caller is responsible for freeing the `result.rates` array.
 */
static XirrCalcResult calculate_xirr_core(XirrWindowContext *ctx) {
    XirrCalcResult result = {NULL, 0};
    double *temp_dates, *temp_values;
    int total_count;

    // Prepare temporary arrays for calculations. This also validates that cash flows are usable.
    if (!prepare_temp_cashflow_data(ctx, &temp_dates, &temp_values, &total_count)) {
        return result; // Return empty result if data is invalid or insufficient.
    }

    // Determine which set of initial guesses to use.
    double default_initial_guesses[] = DEFAULT_INITIAL_GUESSES;
    double *starting_rates = default_initial_guesses;
    int num_starts = sizeof(default_initial_guesses) / sizeof(default_initial_guesses[0]);

    // If the user provided custom starting rates, use them instead of the defaults.
    if (ctx->custom_starting_rates && ctx->num_custom_rates > 0) {
        starting_rates = ctx->custom_starting_rates;
        num_starts = ctx->num_custom_rates;
    }

    // Allocate memory to store the results from each starting guess.
    result.rates = (double *)malloc(num_starts * sizeof(double));
    if (!result.rates) {
        free(temp_values);
        free(temp_dates);
        return result; // Memory allocation failure.
    }
    result.count = num_starts;

    // Run the Newton-Raphson algorithm for each starting guess.
    for (int i = 0; i < num_starts; i++) {
        result.rates[i] = find_root_newton_raphson(starting_rates[i], temp_dates, temp_values, total_count);
    }

    // Clean up the temporary arrays.
    free(temp_values);
    free(temp_dates);
    return result;
}

/**
 * @brief Core engine for `xirr_scan`.
 * @param ctx The window context with all cash flow data.
 * @return A dynamically allocated string of unique, sorted rates. The caller must free this string.
 */
static char *calculate_xirr_scan(XirrWindowContext *ctx) {
    double *temp_dates, *temp_values;
    int total_count;

    // Prepare and validate cash flow data.
    if (!prepare_temp_cashflow_data(ctx, &temp_dates, &temp_values, &total_count)) {
        return NULL;
    }

    // Set scan parameters, using defaults or parsing a user-provided config string.
    double scan_start = DEFAULT_SCAN_START, scan_end = DEFAULT_SCAN_END, scan_step = DEFAULT_SCAN_STEP;
    if (ctx->scan_config_str) {
        parse_scan_config(NULL, ctx->scan_config_str, &scan_start, &scan_end, &scan_step);
    }

    // Dynamically growing array to store the roots (IRRs) we find.
    double *found_rates = NULL;
    int found_count = 0, found_capacity = 0;
    double prev_npv = calculate_npv(scan_start, temp_dates, temp_values, total_count);

    // Scan a range of interest rates to find potential IRR solutions.
    // The strategy is to detect a sign change in the Net Present Value (NPV)
    // between two consecutive steps. A sign change (e.g., from + to -) implies
    // that a root of the NPV function, which is an IRR, exists between those two points.
    for (double r = scan_start + scan_step; r <= scan_end; r += scan_step) {
        double current_npv = calculate_npv(r, temp_dates, temp_values, total_count);
        // A sign change is detected if the product of two consecutive NPVs is negative.
        if (prev_npv * current_npv < 0) {
            // Use Newton-Raphson to precisely locate the root within this small interval.
            // We start the search from the midpoint of the interval where the sign change occurred.
            double root = find_root_newton_raphson(r - (scan_step / 2.0), temp_dates, temp_values, total_count);

            // Check if the found root is a valid, finite number within our defined bounds.
            if (!isnan(root) && !isinf(root) && root > XIRR_MIN_RATE && root < XIRR_MAX_RATE) {
                // Before adding the root, ensure it's not a duplicate of one we've already found.
                int is_duplicate = 0;
                for (int i = 0; i < found_count; i++) {
                    if (fabs(root - found_rates[i]) < XIRR_TOLERANCE) {
                        is_duplicate = 1;
                        break;
                    }
                }
                if (!is_duplicate) {
                    // If the array is full, grow it.
                    if (found_count >= found_capacity) {
                        found_capacity = (found_capacity == 0) ? 10 : found_capacity * 2;
                        double *new_rates = realloc(found_rates, found_capacity * sizeof(double));
                        if (!new_rates) {
                            free(found_rates);
                            free(temp_values);
                            free(temp_dates);
                            return NULL; // Memory allocation failure.
                        }
                        found_rates = new_rates;
                    }
                    // Add the new, unique root to our list.
                    found_rates[found_count++] = root;
                }
            }
        }
        // The current NPV becomes the previous NPV for the next iteration.
        prev_npv = current_npv;
    }

    // Clean up temporary data arrays.
    free(temp_values);
    free(temp_dates);
    // Format the list of found rates into a single string for SQLite.
    char *result_string = format_rates_to_string(found_rates, found_count, 0, 0, 1);
    free(found_rates);
    return result_string;
}

// --- SQLite Callback Functions ---

/**
 * @brief Unified step function for `xirr`, `xirr_all`, `xirr_unique`.
 * @param context The SQLite function context.
 * @param argc The number of arguments passed to the function.
 * @param argv The array of argument values.
 */
static void xirr_step_unified(sqlite3_context *context, int argc, sqlite3_value **argv) {
    if (argc < 2 || argc > 4) {
        sqlite3_result_error(context, "XIRR functions require 2-4 args", -1);
        return;
    }

    // Get or create the aggregate context for this function call.
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, sizeof(XirrWindowContext));
    if (!ctx) {
        sqlite3_result_error_nomem(context);
        return;
    }

    // Initialize the context on the first call.
    if (ctx->data.values == NULL) {
        init_xirr_data(context, &ctx->data);
    }

    // Ignore rows with NULL dates or values.
    if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
        return;
    }

    // Grow the circular buffer if it's full.
    if (ctx->data.count >= ctx->data.capacity) {
        if (grow_xirr_buffer(context, &ctx->data) != SQLITE_OK) {
            return;
        }
    }

    // Parse the date and add the cash flow to the buffer.
    double date_days;
    const char *date_str = (const char *)sqlite3_value_text(argv[0]);
    if (parse_date_string(date_str, &date_days) != SQLITE_OK) {
        sqlite3_result_error(context, "Invalid date format", -1);
        return;
    }
    add_to_circular_buffer(&ctx->data, date_days, sqlite3_value_double(argv[1]));

    // Process optional arguments (current_value, custom_guesses).
    for (int i = 2; i < argc; i++) {
        int arg_type = sqlite3_value_type(argv[i]);
        if (arg_type == SQLITE_NULL)
            continue;

        if (arg_type == SQLITE_FLOAT || arg_type == SQLITE_INTEGER) {
            ctx->current_value = sqlite3_value_double(argv[i]);
        } else if (arg_type == SQLITE_TEXT) {
            if (ctx->custom_starting_rates) {
                free(ctx->custom_starting_rates);
                ctx->custom_starting_rates = NULL;
                ctx->num_custom_rates = 0;
            }
            parse_custom_rates(context, (const char *)sqlite3_value_text(argv[i]), &ctx->custom_starting_rates, &ctx->num_custom_rates);
        }
    }
}

/**
 * @brief Step function for `xirr_scan`.
 * @param context The SQLite function context.
 * @param argc The number of arguments passed to the function.
 * @param argv The array of argument values.
 */
static void xirr_scan_step(sqlite3_context *context, int argc, sqlite3_value **argv) {
    if (argc < 2 || argc > 4) {
        sqlite3_result_error(context, "XIRR_SCAN requires 2-4 args", -1);
        return;
    }

    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, sizeof(XirrWindowContext));
    if (!ctx) {
        sqlite3_result_error_nomem(context);
        return;
    }

    if (ctx->data.values == NULL) {
        init_xirr_data(context, &ctx->data);
    }

    if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
        return;
    }

    if (ctx->data.count >= ctx->data.capacity) {
        if (grow_xirr_buffer(context, &ctx->data) != SQLITE_OK) {
            return;
        }
    }

    double date_days;
    const char *date_str = (const char *)sqlite3_value_text(argv[0]);
    if (parse_date_string(date_str, &date_days) != SQLITE_OK) {
        sqlite3_result_error(context, "Invalid date format", -1);
        return;
    }
    add_to_circular_buffer(&ctx->data, date_days, sqlite3_value_double(argv[1]));

    // Process optional arguments (current_value, scan_config).
    for (int i = 2; i < argc; i++) {
        int arg_type = sqlite3_value_type(argv[i]);
        if (arg_type == SQLITE_NULL)
            continue;

        if (arg_type == SQLITE_FLOAT || arg_type == SQLITE_INTEGER) {
            ctx->current_value = sqlite3_value_double(argv[i]);
        } else if (arg_type == SQLITE_TEXT) {
            if (ctx->scan_config_str) {
                free(ctx->scan_config_str);
                ctx->scan_config_str = NULL;
            }
            ctx->scan_config_str = strdup((const char *)sqlite3_value_text(argv[i]));
        }
    }
}

/**
 * @brief Unified inverse function for all XIRR variants (window functions).
 * @param context The SQLite function context.
 * @param argc The number of arguments passed to the function.
 * @param argv The array of argument values for the row leaving the window.
 */
static void xirr_inverse_unified(sqlite3_context *context, int argc, sqlite3_value **argv) {
    // Get the context without allocating new memory.
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    // If the context and data exist, remove the oldest element from the circular buffer.
    if (ctx && ctx->data.values && ctx->data.count > 0) {
        remove_from_circular_buffer(&ctx->data);
        // Reset current_value if the window becomes empty.
        if (ctx->data.count == 0) {
            ctx->current_value = 0.0;
        }
    }
}

/**
 * @brief Final function for `xirr`. Calculates and returns a single, stable rate.
 * @param context The SQLite function context.
 */
static void xirr_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    XirrCalcResult result = calculate_xirr_core(ctx);
    double final_rate = NAN, boundary_rate = NAN;

    if (result.rates) {
        // Iterate through the results to find the best rate.
        for (int i = 0; i < result.count; i++) {
            if (!isnan(result.rates[i]) && !isinf(result.rates[i])) {
                // The first valid, non-boundary rate is considered the most stable.
                if (result.rates[i] > XIRR_MIN_RATE && result.rates[i] < XIRR_MAX_RATE) {
                    final_rate = result.rates[i];
                    break; // Found a good rate, no need to look further.
                } else if (isnan(boundary_rate)) {
                    // If we only find boundary rates, keep the first one as a fallback.
                    boundary_rate = result.rates[i];
                }
            }
        }
        free(result.rates);
    }

    // If no stable rate was found, use the boundary rate if available.
    if (isnan(final_rate)) {
        final_rate = boundary_rate;
    }
    set_result_double(context, final_rate);
}

/**
 * @brief Value function for `xirr`. For window functions, calculates the result for the current frame.
 * @param context The SQLite function context.
 */
static void xirr_value(sqlite3_context *context) { xirr_final(context); }

/**
 * @brief Final function for `xirr_all`. Returns all found rates as a pipe-separated string.
 * @param context The SQLite function context.
 */
static void xirr_all_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    XirrCalcResult result = calculate_xirr_core(ctx);
    char *result_str = format_rates_to_string(result.rates, result.count, 1, 0, 0);

    if (result_str) {
        // Let SQLite manage the memory of the result string.
        sqlite3_result_text(context, result_str, -1, free);
    } else {
        sqlite3_result_null(context);
    }

    if (result.rates) {
        free(result.rates);
    }
}

/**
 * @brief Value function for `xirr_all`.
 * @param context The SQLite function context.
 */
static void xirr_all_value(sqlite3_context *context) { xirr_all_final(context); }

/**
 * @brief Final function for `xirr_unique`. Returns unique, non-boundary rates as a string.
 * @param context The SQLite function context.
 */
static void xirr_unique_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    XirrCalcResult result = calculate_xirr_core(ctx);
    char *result_str = format_rates_to_string(result.rates, result.count, 0, 1, 0);

    if (result_str) {
        sqlite3_result_text(context, result_str, -1, free);
    } else {
        sqlite3_result_null(context);
    }

    if (result.rates) {
        free(result.rates);
    }
}

/**
 * @brief Value function for `xirr_unique`.
 * @param context The SQLite function context.
 */
static void xirr_unique_value(sqlite3_context *context) { xirr_unique_final(context); }

/**
 * @brief Final function for `xirr_scan`. Returns all found rates from scanning, sorted.
 * @param context The SQLite function context.
 */
static void xirr_scan_final(sqlite3_context *context) {
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    char *result_str = calculate_xirr_scan(ctx);

    if (result_str) {
        sqlite3_result_text(context, result_str, -1, free);
    } else {
        sqlite3_result_null(context);
    }
}

/**
 * @brief Value function for `xirr_scan`.
 * @param context The SQLite function context.
 */
static void xirr_scan_value(sqlite3_context *context) { xirr_scan_final(context); }

/**
 * @brief Destructor for the aggregate context.
 *
 * This function is registered with SQLite and is guaranteed to be called,
 * even if the query is aborted or encounters an error. It ensures that all
 * dynamically allocated memory within the context is freed, preventing memory leaks.
 * @param pAggregate The aggregate context to be destroyed.
 */
static void xirr_destroy(void *pAggregate) {
    XirrWindowContext *ctx = (XirrWindowContext *)pAggregate;
    if (ctx) {
        if (ctx->data.dates)
            free(ctx->data.dates);
        if (ctx->data.values)
            free(ctx->data.values);
        if (ctx->custom_starting_rates)
            free(ctx->custom_starting_rates);
        if (ctx->scan_config_str)
            free(ctx->scan_config_str);
    }
}

/**
 * @brief Main entry point for the SQLite extension. Registers all XIRR functions.
 * @param db The database connection handle.
 * @param pzErrMsg A pointer to an error message string.
 * @param pApi A pointer to the SQLite API routines.
 * @return SQLITE_OK on success, or an error code on failure.
 */
int sqlite3_extension_init(sqlite3 *db, char **pzErrMsg, const sqlite3_api_routines *pApi) {
    int rc = SQLITE_OK;
    SQLITE_EXTENSION_INIT2(pApi);
    int flags = SQLITE_UTF8 | SQLITE_DETERMINISTIC | SQLITE_INNOCUOUS;

    // Register xirr and its alias XIRR
    rc = sqlite3_create_window_function(db, "xirr", -1, flags, 0, xirr_step_unified, xirr_final, xirr_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;
    rc = sqlite3_create_window_function(db, "XIRR", -1, flags, 0, xirr_step_unified, xirr_final, xirr_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;

    // Register xirr_all and its alias XIRR_ALL
    rc = sqlite3_create_window_function(db, "xirr_all", -1, flags, 0, xirr_step_unified, xirr_all_final, xirr_all_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;
    rc = sqlite3_create_window_function(db, "XIRR_ALL", -1, flags, 0, xirr_step_unified, xirr_all_final, xirr_all_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;

    // Register xirr_unique and its alias XIRR_UNIQUE
    rc = sqlite3_create_window_function(db, "xirr_unique", -1, flags, 0, xirr_step_unified, xirr_unique_final, xirr_unique_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;
    rc = sqlite3_create_window_function(db, "XIRR_UNIQUE", -1, flags, 0, xirr_step_unified, xirr_unique_final, xirr_unique_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;

    // Register xirr_scan and its alias XIRR_SCAN
    rc = sqlite3_create_window_function(db, "xirr_scan", -1, flags, 0, xirr_scan_step, xirr_scan_final, xirr_scan_value, xirr_inverse_unified, xirr_destroy);
    if (rc != SQLITE_OK)
        return rc;
    rc = sqlite3_create_window_function(db, "XIRR_SCAN", -1, flags, 0, xirr_scan_step, xirr_scan_final, xirr_scan_value, xirr_inverse_unified, xirr_destroy);

    return rc;
}

// --- Helper Function Implementations ---

/**
 * @brief Prepares temporary cash flow data arrays for calculation.
 * @param ctx The window context.
 * @param out_dates Pointer to store the allocated dates array.
 * @param out_values Pointer to store the allocated values array.
 * @param out_count Pointer to store the total count of cash flows.
 * @return 1 on success, 0 on failure (e.g., memory error, invalid data).
 */
static int prepare_temp_cashflow_data(XirrWindowContext *ctx, double **out_dates, double **out_values, int *out_count) {
    if (!ctx || !ctx->data.values || ctx->data.count < 1) {
        return 0;
    }

    // Total count includes all historical values plus the current one.
    int total_count = ctx->data.count + 1;
    double *temp_dates = (double *)malloc(total_count * sizeof(double));
    double *temp_values = (double *)malloc(total_count * sizeof(double));

    if (!temp_dates || !temp_values) {
        if (temp_dates)
            free(temp_dates);
        if (temp_values)
            free(temp_values);
        return 0;
    }

    // Copy data from the circular buffer to the temporary flat arrays.
    for (int i = 0; i < ctx->data.count; i++) {
        get_circular_data(&ctx->data, i, &temp_dates[i], &temp_values[i]);
    }
    // Add the final cash flow value for the current window frame.
    temp_values[ctx->data.count] = ctx->current_value;
    double last_date, last_val;
    get_circular_data(&ctx->data, ctx->data.count - 1, &last_date, &last_val);
    temp_dates[ctx->data.count] = last_date;

    // The XIRR algorithm requires at least one positive and one negative cash flow.
    int has_positive = 0, has_negative = 0;
    for (int i = 0; i < total_count; i++) {
        if (temp_values[i] > 0)
            has_positive = 1;
        if (temp_values[i] < 0)
            has_negative = 1;
    }

    // If this condition is not met, a meaningful IRR cannot be calculated.
    if (!has_positive || !has_negative) {
        free(temp_values);
        free(temp_dates);
        return 0;
    }

    *out_dates = temp_dates;
    *out_values = temp_values;
    *out_count = total_count;
    return 1;
}

/**
 * @brief Checks if a year is a leap year.
 * @param year The year to check.
 * @return 1 if leap, 0 otherwise.
 */
static int is_leap_year(int year) { return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0); }

/**
 * @brief Returns the number of days in a given month of a year.
 * @param year The year.
 * @param month The month (1-12).
 * @return The number of days.
 */
static int days_in_month(int year, int month) {
    int d[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return (month == 2 && is_leap_year(year)) ? 29 : d[month - 1];
}

/**
 * @brief Converts a date to a number of days since MIN_DATE_YEAR.
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
 * @brief Parses a date string (e.g., "YYYY-MM-DD") into days.
 * @param date_str The date string.
 * @param out_days Pointer to store the result.
 * @return SQLITE_OK on success.
 */
static int parse_date_string(const char *date_str, double *out_days) {
    if (!date_str)
        return SQLITE_MISMATCH;
    int y, m, d;
    // Try parsing YYYY-MM-DD format
    if (sscanf(date_str, "%d-%d-%d", &y, &m, &d) == 3) {
        if (y >= MIN_DATE_YEAR && y <= 9999 && m >= 1 && m <= 12 && d >= 1 && d <= days_in_month(y, m)) {
            *out_days = date_to_days_from_components(y, m, d);
            return SQLITE_OK;
        }
    }

    return SQLITE_MISMATCH;
}

/**
 * @brief Parses a pipe-delimited string of initial guess rates.
 * @param context SQLite context.
 * @param rates_str The string to parse.
 * @param out_rates Pointer to the allocated array of rates.
 * @param out_count Pointer to the number of rates.
 * @return SQLITE_OK on success.
 */
static int parse_custom_rates(sqlite3_context *context, const char *rates_str, double **out_rates, int *out_count) {
    if (!rates_str || *rates_str == '\0') {
        return SQLITE_MISMATCH;
    }

    // First pass: count the number of rates to allocate the correct amount of memory.
    int count = 1;
    const char *p = rates_str;
    while (*p) {
        if (*p++ == '|') {
            count++;
        }
    }

    double *rates = (double *)malloc(count * sizeof(double));
    if (!rates) {
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }

    // Second pass: parse the rates using strtod for safe conversion.
    const char *start = rates_str;
    int i = 0;
    while (start && *start && i < count) {
        char *endptr;
        rates[i] = strtod(start, &endptr);

        // Check if a valid number was parsed.
        if (endptr == start) {
            free(rates);
            sqlite3_result_error(context, "Invalid number in guesses string", -1);
            return SQLITE_MISMATCH;
        }

        // Check for invalid characters after the number.
        if (*endptr != '\0' && *endptr != '|') {
            free(rates);
            sqlite3_result_error(context, "Invalid characters in guesses string", -1);
            return SQLITE_MISMATCH;
        }

        i++;
        // Move the start pointer to the beginning of the next number.
        start = (*endptr == '|') ? (endptr + 1) : NULL;
    }

    *out_rates = rates;
    *out_count = i;
    return SQLITE_OK;
}

/**
 * @brief Parses a pipe-delimited string for scan configuration (`start|end|step`).
 * @param context SQLite context.
 * @param config_str The string to parse.
 * @param start Pointer to the start rate.
 * @param end Pointer to the end rate.
 * @param step Pointer to the step size.
 * @return SQLITE_OK on success.
 */
static int parse_scan_config(sqlite3_context *context, const char *config_str, double *start, double *end, double *step) {
    if (!config_str || *config_str == '\0') {
        return SQLITE_OK; // Nothing to parse, use defaults.
    }

    const char *p = config_str;
    char *endptr;
    double val;

    // Parse start value.
    val = strtod(p, &endptr);
    if (endptr == p)
        return SQLITE_OK; // Empty or invalid string, use defaults.
    *start = val;
    p = endptr;

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
    if (*p != '\0') {
        if (context) {
            sqlite3_result_error(context, "Invalid characters in scan config string", -1);
        }
        return SQLITE_MISMATCH;
    }

    return SQLITE_OK;
}

/**
 * @brief Calculates the Net Present Value (NPV).
 * @param rate The discount rate.
 * @param dates Array of dates.
 * @param values Array of cash flows.
 * @param count Number of cash flows.
 * @return The calculated NPV.
 */
static double calculate_npv(double rate, double *dates, double *values, int count) {
    if (count < 2)
        return NAN;

    double npv = 0.0;
    double base_date = dates[0];
    for (int i = 0; i < count; i++) {
        // Calculate the time difference in years from the first cash flow.
        double years = (dates[i] - base_date) / DAYS_PER_YEAR;
        // Calculate the discount factor.
        double factor = pow(1.0 + rate, years);
        if (isinf(factor) || isnan(factor))
            return NAN;
        // Add the discounted cash flow to the total NPV.
        npv += values[i] / factor;
    }
    return npv;
}

/**
 * @brief Calculates the derivative of the NPV function.
 * @param rate The discount rate.
 * @param dates Array of dates.
 * @param values Array of cash flows.
 * @param count Number of cash flows.
 * @return The calculated derivative.
 */
static double calculate_npv_derivative(double rate, double *dates, double *values, int count) {
    if (count < 2)
        return NAN;

    double derivative = 0.0;
    double base_date = dates[0];
    for (int i = 0; i < count; i++) {
        double years = (dates[i] - base_date) / DAYS_PER_YEAR;
        double factor = pow(1.0 + rate, years + 1.0);
        if (isinf(factor) || isnan(factor))
            return NAN;
        derivative -= years * values[i] / factor;
    }
    return derivative;
}

/**
 * @brief Retrieves a date/value pair from the circular buffer.
 * @param data The XirrData struct.
 * @param logical_index The logical index to retrieve.
 * @param date Pointer to store the date.
 * @param value Pointer to store the value.
 */
static void get_circular_data(XirrData *data, int logical_index, double *date, double *value) {
    // Map the logical index to the physical index in the circular buffer array.
    int phys_idx = (data->head + logical_index) % data->capacity;
    *date = data->dates[phys_idx];
    *value = data->values[phys_idx];
}

/**
 * @brief Adds a new cash flow to the circular buffer.
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
 * @brief Removes the oldest cash flow from the circular buffer.
 * @param data The XirrData struct.
 */
static void remove_from_circular_buffer(XirrData *data) {
    if (data->count == 0)
        return;
    // Move the head forward, effectively "removing" the oldest element.
    data->head = (data->head + 1) % data->capacity;
    data->count--;
}

/**
 * @brief Initializes the XirrData structure.
 * @param context SQLite context.
 * @param data The XirrData struct to initialize.
 * @return SQLITE_OK on success.
 */
static int init_xirr_data(sqlite3_context *context, XirrData *data) {
    data->capacity = INITIAL_CAPACITY;
    data->count = 0;
    data->head = 0;
    data->tail = 0;
    data->dates = (double *)malloc(data->capacity * sizeof(double));
    data->values = (double *)malloc(data->capacity * sizeof(double));
    if (!data->dates || !data->values) {
        if (data->dates)
            free(data->dates);
        if (data->values)
            free(data->values);
        sqlite3_result_error_nomem(context);
        return SQLITE_NOMEM;
    }
    // Initialize all other context fields to zero/NULL.
    XirrWindowContext *ctx = (XirrWindowContext *)sqlite3_aggregate_context(context, 0);
    if (ctx) {
        ctx->current_value = 0.0;
        ctx->custom_starting_rates = NULL;
        ctx->num_custom_rates = 0;
        ctx->scan_config_str = NULL;
    }
    return SQLITE_OK;
}

/**
 * @brief Grows the circular buffer.
 * @param context SQLite context.
 * @param data The XirrData struct to grow.
 * @return SQLITE_OK on success.
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

    // Copy existing data from the old circular buffer to the new linear arrays.
    for (int i = 0; i < data->count; i++) {
        get_circular_data(data, i, &new_dates[i], &new_values[i]);
    }

    free(data->dates);
    free(data->values);

    // Update the data pointers and reset head/tail for the new, larger linear layout.
    data->dates = new_dates;
    data->values = new_values;
    data->capacity = new_capacity;
    data->head = 0;
    data->tail = data->count;
    return SQLITE_OK;
}

/**
 * @brief Sets the SQLite result to a double, handling NAN/INF.
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
 * @return Comparison result.
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
 * @param rates Array of rates.
 * @param count Number of rates.
 * @param include_empty If true, include non-converged rates as empty strings.
 * @param filter_unique If true, filter for unique, non-boundary rates.
 * @param sort_desc If true, sort results descending.
 * @return A dynamically allocated string.
 */
static char *format_rates_to_string(double *rates, int count, int include_empty, int filter_unique, int sort_desc) {
    if (!rates || count == 0)
        return NULL;

    if (sort_desc) {
        qsort(rates, count, sizeof(double), compare_doubles_desc);
    }

    // Create a temporary array to hold the rates that will be included in the final string.
    double *final_rates = (double *)malloc(count * sizeof(double));
    if (!final_rates)
        return NULL;
    int final_count = 0;

    // First pass: filter the rates according to the function's parameters.
    for (int i = 0; i < count; i++) {
        if (isnan(rates[i]) || isinf(rates[i])) {
            if (include_empty) {
                final_rates[final_count++] = NAN;
            }
            continue;
        }

        if (filter_unique) {
            if (rates[i] <= XIRR_MIN_RATE || rates[i] >= XIRR_MAX_RATE) {
                continue; // Exclude boundary rates.
            }
            int is_duplicate = 0;
            for (int j = 0; j < final_count; j++) {
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
        free(final_rates);
        return NULL;
    }

    // Second pass: calculate the exact size needed for the result string.
    size_t total_len = 0;
    for (int i = 0; i < final_count; i++) {
        if (!isnan(final_rates[i])) {
            char buffer[64];
            total_len += snprintf(buffer, sizeof(buffer), "%.15g", final_rates[i]);
        }
    }
    total_len += (final_count > 0) ? (final_count - 1) : 0; // Add space for separators.
    total_len += 1;                                         // Add space for the null terminator.

    // Third pass: allocate memory and build the string.
    char *result_string = (char *)malloc(total_len);
    if (!result_string) {
        free(final_rates);
        return NULL;
    }

    char *ptr = result_string;
    size_t remaining_space = total_len;
    for (int i = 0; i < final_count; i++) {
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
    *ptr = '\0'; // Null-terminate the string.

    free(final_rates);
    return result_string;
}

/**
 * @brief Runs the Newton-Raphson algorithm to find a precise root of the NPV function.
 *
 * The method iteratively improves the guess for the rate using the formula:
 *   new_rate = rate - NPV(rate) / NPV'(rate)
 * where NPV' is the derivative of the NPV function. The process stops when the
 * change in the rate is below the tolerance threshold, the NPV is close to zero,
 * or the maximum number of iterations is reached.
 *
 * @param start_rate The initial guess.
 * @param dates Array of dates.
 * @param values Array of cash flows.
 * @param count Number of cash flows.
 * @return The converged rate, or NAN if it fails.
 */
static double find_root_newton_raphson(double start_rate, double *dates, double *values, int count) {
    double rate = start_rate;
    for (int i = 0; i < XIRR_MAX_ITERATIONS; i++) {
        double npv = calculate_npv(rate, dates, values, count);
        double npv_derivative = calculate_npv_derivative(rate, dates, values, count);

        // Stop if the calculation results in non-finite numbers or if the derivative is too close to zero.
        if (isnan(npv) || isinf(npv) || isnan(npv_derivative) || isinf(npv_derivative) || fabs(npv_derivative) < XIRR_TOLERANCE) {
            break;
        }
        // If the NPV is very close to zero, we have found the root.
        if (fabs(npv) < XIRR_TOLERANCE) {
            return rate;
        }

        // Calculate the next guess for the rate.
        double new_rate = rate - npv / npv_derivative;
        // Clamp the new rate to our predefined min/max bounds.
        if (new_rate < XIRR_MIN_RATE)
            new_rate = XIRR_MIN_RATE;
        if (new_rate > XIRR_MAX_RATE)
            new_rate = XIRR_MAX_RATE;

        // If the change in the rate is negligible, we have converged.
        if (fabs(new_rate - rate) < XIRR_TOLERANCE) {
            return new_rate;
        }
        rate = new_rate;
    }
    return NAN; // Return NAN if the algorithm did not converge.
}