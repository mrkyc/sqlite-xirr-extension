# XIRR Extension for SQLite

This code implements a suite of `XIRR` (Extended Internal Rate of Return) functions as a SQLite extension. It allows you to calculate the internal rate of return for a series of irregular cash flows, providing multiple ways to handle complex scenarios where more than one solution to the XIRR equation may exist.

This extension provides four main functions:
-   `xirr()`: Calculates a single, stable XIRR value. It is designed to find the most likely rate of return and is suitable for most standard financial calculations.
-   `xirr_all()`: Returns all potential XIRR solutions found based on a set of initial guesses.
-   `xirr_unique()`: Filters the results from `xirr_all` to return only the unique solutions.
-   `xirr_scan()`: Performs a systematic scan across a wide range of rates to find all possible solutions, which are then returned sorted in descending order.

All functions are available as both aggregate and window functions and are registered in both lowercase (`xirr`) and uppercase (`XIRR`).

## Table of Contents

- [How It Works](#how-it-works)
- [Available Functions](#available-functions)
- [Compilation and Loading](#compilation-and-loading)
- [Usage](#usage)
  - [Syntax](#syntax)
  - [Arguments](#arguments)
- [Examples](#examples)
  - [SQL Example (Aggregate)](#sql-example-aggregate)
  - [SQL Example (Window Function)](#sql-example-window-function)
  - [Python Example](#python-example)
- [Limitations and Error Handling](#limitations-and-error-handling)

## How It Works

The core of the extension uses two primary methods to find the discount rate(s) at which the Net Present Value (NPV) of a series of cash flows equals zero:

1.  **Newton-Raphson Method (`xirr`, `xirr_all`, `xirr_unique`):** This is a fast, iterative root-finding algorithm. The functions start from a list of initial guesses and converge on a solution. This method is efficient but may not find all possible solutions if the initial guesses are not well-chosen.
2.  **Scanning with Newton-Raphson Refinement (`xirr_scan`):** This method systematically scans a wide range of rates (e.g., from -99% to 1000%) with a small step. When it detects a sign change in the NPV between two steps (which implies a root exists in that interval), it uses the Newton-Raphson method to pinpoint its exact value. This approach is more robust for finding all solutions but is computationally more intensive.

## Available Functions

### `xirr(date, cash_flow, [arg3], [arg4], [arg5])`
-   **Returns:** A single floating-point number (`DOUBLE`).
-   **Description:** Searches for a stable XIRR using the Newton-Raphson method. It returns the first valid solution found from its list of initial guesses. Ideal for standard calculations.

### `xirr_all(date, cash_flow, [arg3], [arg4], [arg5])`
-   **Returns:** A pipe-separated string (`TEXT`).
-   **Description:** For each provided initial guess, this function runs the calculation and returns the corresponding result. This includes all valid rates and duplicates. If a calculation for a specific initial guess does not converge to a valid solution (which includes rates that fall exactly on the boundaries), its place in the output string will be empty, preserving the order (e.g., `0.25||-0.1` indicates the second guess failed).

### `xirr_unique(date, cash_flow, [arg3], [arg4], [arg5])`
-   **Returns:** A pipe-separated string (`TEXT`).
-   **Description:** Works like `xirr_all` but filters the output to include only unique solutions.

### `xirr_scan(date, cash_flow, [arg3], [arg4], [arg5])`
-   **Returns:** A pipe-separated string (`TEXT`).
-   **Description:** Systematically scans a range of rates to find all possible solutions. The results are guaranteed to be unique and are returned sorted in descending order. This is the most reliable function for complex cash flows with multiple roots.

## Compilation and Loading

To use this extension, you first need to compile it into a shared library.

### Prerequisites

Before compiling, you must have the SQLite development libraries installed, which provide the necessary header files (`sqlite3.h` and `sqlite3ext.h`).

-   **On Debian/Ubuntu:** `sudo apt-get install sqlite3-dev`
-   **On Fedora/CentOS:** `sudo dnf install sqlite-devel`
-   **On macOS (with Homebrew):** `brew install sqlite`
-   **On Windows:** Ensure your compiler (like MinGW-w64) has access to the SQLite source code or pre-compiled headers.

### Compilation Instructions

Compile `sqlite-xirr-extension.c` into a shared library.

- **Linux:** `gcc -shared -fPIC -o sqlite-xirr-extension.so sqlite-xirr-extension.c -lm`
- **macOS:** `gcc -shared -fPIC -I$(brew --prefix sqlite)/include -undefined dynamic_lookup -o sqlite-xirr-extension.dylib sqlite-xirr-extension.c -lm`
- **Windows:** `gcc -shared -o sqlite-xirr-extension.dll sqlite-xirr-extension.c -lm`

### Loading the Extension

Once compiled, you can load the extension in your SQLite session:

```sql
-- On Linux
.load ./sqlite-xirr-extension.so

-- On macOS
.load ./sqlite-xirr-extension.dylib

-- On Windows
.load ./sqlite-xirr-extension.dll
```

## Usage

### Syntax

All functions accept between 2 and 5 arguments.

```sql
FUNCTION_NAME(date, cash_flow, [arg3], [arg4], [arg5])
```

### Arguments

1.  `date` (text, **required**): The date of the cash flow in `YYYY-MM-DD` format.
2.  `cash_flow` (numeric, **required**): The value of the cash flow.
3.  **Optional Arguments**: The function can take up to three optional arguments. It intelligently determines their meaning based on type and order:
    *   `current_value` (numeric): The total market value of the investment at the time of the last transaction. This is treated as a final, positive cash flow occurring on the same date as the last transaction in the set. Defaults to `0`.
    *   `start_date` (text): A date string in `YYYY-MM-DD` format. This argument acts as a **calculation trigger, not a data filter**. The function will still consider **all cash flows within the current aggregate or window frame**. However, it will only perform the calculation and return a value if the date of the *last transaction* in that frame is on or after the provided `start_date`. For all rows where this condition is not met, the function will return `NULL`. This is useful for calculating rolling returns that "activate" only after a certain point in time.
    *   `config_string` (text): A configuration string whose meaning depends on the function:
        -   For `xirr`, `xirr_all`, `xirr_unique`: This is the `initial_guesses` string, a pipe-separated list of starting rates (e.g., `'0.1|0.25|-0.05'`). The default is `'0.1|0.05|0.2|-0.1|0.01|0.5'`.
        -   For `xirr_scan`: This is the `scan_config` string, with positional values for `start|end|step` (e.g., `'-0.5|5.0|0.005'`). Defaults are `-0.99` (start), `10.0` (end), and `0.01` (step).

    **Argument Parsing Logic:** The function parses the optional arguments in the order they are given. A numeric argument is always treated as `current_value`. A text argument is first tested to see if it's a valid date; if so, it's treated as `start_date`. If it's not a valid date, it's treated as `config_string`.

## Examples

### SQL Example (Aggregate)

Let's assume we have a `transactions` table. The following examples show how to use each function with different argument combinations.

```sql
-- 2 args: Basic calculation
SELECT XIRR(transaction_date, amount) FROM transactions;

-- 3 args: With current market value
SELECT XIRR(transaction_date, amount, portfolio_value) FROM transactions;

-- 3 args: With custom initial guesses
SELECT XIRR_ALL(transaction_date, amount, '0.15|0.2') FROM transactions;

-- 4 args: With current value and custom guesses
SELECT XIRR_UNIQUE(transaction_date, amount, portfolio_value, '0.15|0.2') FROM transactions;

-- 5 args: With current value, start date, and custom guesses
SELECT XIRR(transaction_date, amount, portfolio_value, '2023-06-01', '0.1|0.2')
FROM transactions;
```

### SQL Example (Window Function)

The functions can also be used over a sliding window, which is useful for calculating rolling returns.

```sql
-- Calculate a rolling 6-month XIRR for each transaction,
-- but only if the transaction is after 2023-03-01
SELECT
  transaction_date,
  amount,
  XIRR(transaction_date, amount, '2023-03-01') OVER (
    ORDER BY transaction_date
    ROWS BETWEEN 5 PRECEDING AND CURRENT ROW
  ) AS rolling_6mo_xirr
FROM transactions;
```

### Python Example

This example demonstrates loading the extension and using the functions.

```python
import sqlite3

# Adjust for your OS and compiled filename
EXTENSION_PATH = './sqlite-xirr-extension.so' 

con = sqlite3.connect(':memory:')
con.enable_load_extension(True)
con.load_extension(EXTENSION_PATH)

print("Extension loaded successfully.")

cur = con.cursor()
cur.execute("CREATE TABLE transactions(portfolio_id, transaction_date, amount, portfolio_value)")
transactions = [
    (1, '2023-01-01', -1000, 1000),
    (1, '2023-03-15', -500,  1600),
    (1, '2023-06-01', 200,   1900),
    (1, '2023-08-20', 100,   2100),
    (2, '2023-01-01', -2000, 2000),
    (2, '2023-04-01', 500,   2500),
    (2, '2023-07-01', 150,   3000)
]
cur.executemany("INSERT INTO transactions VALUES(?, ?, ?, ?)", transactions)

query = """
SELECT
  XIRR(transaction_date, amount, portfolio_value) AS xirr_rate,
  XIRR_ALL(transaction_date, amount, portfolio_value, '0.1|0.5') AS all_rates,
  XIRR_UNIQUE(transaction_date, amount, portfolio_value, '0.1|0.5') AS unique_rates,
  XIRR_SCAN(transaction_date, amount, portfolio_value, '-0.5|2.0|0.01') AS scan_rates
FROM transactions
GROUP BY portfolio_id;
"""

print("--- XIRR Results ---")
for row in cur.execute(query):
    print(f"XIRR: {row[0]}
All Rates: {row[1]}
Unique Rates: {row[2]}
Scan Rates: {row[3]}")

con.close()
```

## Limitations and Error Handling

-   **Date Range:** Supports dates from **1900-01-01** to **9999-12-31**.
-   **Prerequisite:** Requires at least one positive and one negative cash flow. If this condition is not met, the functions will return `NULL`.
-   **Rate Boundaries:** All successfully calculated rates are guaranteed to be within the predefined bounds: greater than -0.9999 (-99.99%) and less than 10.0 (1000%).
-   **Convergence:** The Newton-Raphson method (`xirr`, `xirr_all`, `xirr_unique`) may not always converge. In such cases, `xirr` will return `NULL`, and `xirr_all` will have an empty section in its output string. `xirr_scan` is more reliable for finding all existing roots within its scan range.
-   **Performance:** `xirr_scan` is significantly more resource-intensive than the other functions due to its comprehensive search.
-   **Internal Error Handling (C Code):**
    -   **Invalid Arguments:** The C code performs rigorous validation of input arguments. It checks for the correct data types (e.g., `TEXT` for dates, `NUMERIC` for cash flows) and valid formats (e.g., `YYYY-MM-DD` for dates, pipe-separated numbers for custom guesses/scan config). If an argument is invalid, `sqlite3_result_error` is used to return a descriptive error message to SQLite, halting query execution.
    -   **Memory Allocation Failures:** The extension dynamically allocates memory for cash flow data and other internal structures. If a memory allocation fails (e.g., `malloc` or `realloc` returns `NULL`), `sqlite3_result_error_nomem` is called to signal an out-of-memory condition to SQLite.
    -   **Calculation Edge Cases & Non-Convergence:**
        -   **Insufficient Data:** Functions return `NULL` if there are no cash flows or if the cash flows do not contain both positive and negative values (a prerequisite for XIRR calculation).
        -   **Non-Finite Results:** If intermediate or final calculations result in `NaN` (Not a Number) or `Infinity` (e.g., due to division by zero or extreme values), `sqlite3_result_null` is used to return `NULL`.
        -   **Newton-Raphson Convergence:** The Newton-Raphson method used by `xirr`, `xirr_all`, and `xirr_unique` may not always converge to a solution within the maximum number of iterations. In such cases, `NAN` is returned internally, which then translates to `NULL` in SQLite results for `xirr`, or an empty section in the pipe-separated string for `xirr_all` and `xirr_unique`.