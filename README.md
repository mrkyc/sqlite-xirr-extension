# XIRR Extension for SQLite

This code implements a suite of `XIRR` (Extended Internal Rate of Return) functions as a SQLite extension. It allows you to calculate the internal rate of return for a series of irregular cash flows, providing multiple ways to handle complex scenarios where more than one solution to the XIRR equation may exist.

This extension provides four main functions:
-   `xirr()`: Calculates a single, stable XIRR value. It is designed to find the most likely rate of return and is suitable for most standard financial calculations.
-   `xirr_all()`: Returns all potential XIRR solutions found based on a set of initial guesses.
-   `xirr_unique()`: Filters the results from `xirr_all` to return only the unique, non-boundary solutions.
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
  - [SQL Example](#sql-example)
  - [Python Example](#python-example)
- [Limitations](#limitations)

## How It Works

The core of the extension uses two primary methods to find the discount rate(s) at which the Net Present Value (NPV) of a series of cash flows equals zero:

1.  **Newton-Raphson Method (`xirr`, `xirr_all`, `xirr_unique`):** This is a fast, iterative root-finding algorithm. The functions start from a list of initial guesses and converge on a solution. This method is efficient but may not find all possible solutions if the initial guesses are not well-chosen.
2.  **Bisection and Scanning (`xirr_scan`):** This method systematically scans a wide range of rates (e.g., from -99% to 1000%) with a small step. When it detects a sign change in the NPV between two steps, it knows a root exists in that interval and then uses the Newton-Raphson method to pinpoint its exact value. This approach is more robust for finding all solutions but is computationally more intensive.

## Available Functions

### `xirr(date, cash_flow, [current_value], [initial_guesses])`
-   **Returns:** A single floating-point number (`DOUBLE`).
-   **Description:** Searches for a stable XIRR using the Newton-Raphson method. It returns the first valid, non-boundary solution found from its list of initial guesses. Ideal for standard calculations.

### `xirr_all(date, cash_flow, [current_value], [initial_guesses])`
-   **Returns:** A pipe-separated string (`TEXT`).
-   **Description:** Tests every initial guess and returns all resulting rates, including duplicates and boundary solutions. An empty string is returned for guesses that do not converge.

### `xirr_unique(date, cash_flow, [current_value], [initial_guesses])`
-   **Returns:** A pipe-separated string (`TEXT`).
-   **Description:** Works like `xirr_all` but filters the output to include only unique, non-boundary solutions.

### `xirr_scan(date, cash_flow, [current_value], [scan_config])`
-   **Returns:** A pipe-separated string (`TEXT`).
-   **Description:** Systematically scans a range of rates to find all possible solutions. The results are guaranteed to be unique and are returned sorted in descending order. This is the most reliable function for complex cash flows with multiple roots.

## Compilation and Loading

To use this extension, compile `sqlite-xirr-extension.c` into a shared library.

- **Prerequisites:** SQLite development libraries (`sqlite3.h` and `sqlite3ext.h`).
- **Linux:** `gcc -shared -fPIC -o sqlite-xirr-extension.so sqlite-xirr-extension.c -lm`
- **macOS:** `gcc -shared -fPIC -I$(brew --prefix sqlite)/include -undefined dynamic_lookup -o sqlite-xirr-extension.dylib sqlite-xirr-extension.c -lm`
- **Windows:** `gcc -shared -o sqlite-xirr-extension.dll sqlite-xirr-extension.c -lm`

Load the extension in SQLite:
```sql
.load ./sqlite-xirr-extension.so
```

## Usage

### Syntax

All functions accept 2, 3, or 4 arguments.

```sql
FUNCTION_NAME(date, cash_flow, [current_value], [config_string])
```

### Arguments

1.  `date` (text, **required**): The date of the cash flow in `YYYY-MM-DD` format.
2.  `cash_flow` (numeric, **required**): The value of the cash flow.
3.  `current_value` (numeric, *optional*): The total market value of the investment at the time of the last transaction. Defaults to `0`.
4.  `config_string` (text, *optional*): A configuration string. Its meaning depends on the function:
    -   For `xirr`, `xirr_all`, `xirr_unique`: This is the `initial_guesses` string, a pipe-separated list of starting rates (e.g., `'0.1|0.25|-0.05'`). The default is `'0.1|0.05|0.2|-0.1|0.01|0.5'`.
    -   For `xirr_scan`: This is the `scan_config` string, with positional values for `start|end|step` (e.g., `'-0.5|5.0|0.005'`). Defaults are `-0.99` (start), `10.0` (end), and `0.01` (step).

The function intelligently distinguishes between `current_value` and the `config_string` based on their data type (numeric vs. text).

## Examples

### SQL Example

Let's assume we have a `transactions` table. The following examples show how to use each function with different argument combinations.

```sql
-- 2 args: Basic calculation
SELECT XIRR(transaction_date, amount) FROM transactions;
SELECT XIRR_ALL(transaction_date, amount) FROM transactions;
SELECT XIRR(transaction_date, amount) FROM transactions;
SELECT XIRR_SCAN(transaction_date, amount) FROM transactions;

-- 3 args: With current market value
SELECT XIRR(transaction_date, amount, portfolio_value) FROM transactions;
SELECT XIRR_ALL(transaction_date, amount, portfolio_value) FROM transactions;
SELECT XIRR_UNIQUE(transaction_date, amount, portfolio_value) FROM transactions;
SELECT XIRR_SCAN(transaction_date, amount, portfolio_value) FROM transactions;

-- 3 args: With custom initial guesses
SELECT XIRR(transaction_date, amount, '0.15|0.2') FROM transactions;
SELECT XIRR_ALL(transaction_date, amount, '0.15|0.2') FROM transactions;
SELECT XIRR_UNIQUE(transaction_date, amount, '0.15|0.2') FROM transactions;

-- 3 args: Scan with custom settings (start=-0.5, end=5.0, step=0.005)
SELECT XIRR_SCAN(transaction_date, amount, '-0.5|5.0|0.005') FROM transactions;

-- 4 args: With current value and custom guesses
SELECT XIRR(transaction_date, amount, portfolio_value, '0.15|0.2') FROM transactions;
SELECT XIRR_ALL(transaction_date, amount, portfolio_value, '0.15|0.2') FROM transactions;
SELECT XIRR_UNIQUE(transaction_date, amount, portfolio_value, '0.15|0.2') FROM transactions;

-- 4 args: Scan with current value and custom settings
SELECT XIRR_SCAN(transaction_date, amount, portfolio_value, '-0.5|5.0|0.005') FROM transactions;
```

### Python Example

This example demonstrates loading the extension and using all four functions.

```python
import sqlite3

EXTENSION_PATH = './sqlite-xirr-extension.so' # Adjust for your OS

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

print("\n--- XIRR Results ---")
for row in cur.execute(query):
    print(f"XIRR: {row[0]}\nAll Rates: {row[1]}\nUnique Rates: {row[2]}\nScan Rates: {row[3]}")

con.close()
```

## Limitations

-   **Date Range:** Supports dates from **1900-01-01** to **9999-12-31**.
-   **Prerequisite:** Requires at least one positive and one negative cash flow.
-   **Convergence:** The Newton-Raphson method (`xirr`, `xirr_all`, `xirr_unique`) may not always converge. `xirr_scan` is more reliable for finding all existing roots within its scan range.
-   **Performance:** `xirr_scan` is significantly more resource-intensive than the other functions due to its comprehensive search.
