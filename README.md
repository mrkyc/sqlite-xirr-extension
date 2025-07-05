# XIRR Extension for SQLite

This code implements an `XIRR` (Extended Internal Rate of Return) window function as an extension for the SQLite database. It allows you to calculate the internal rate of return for a series of irregular cash flows.

## Table of Contents

- [How It Works](#how-it-works)
- [Compilation and Loading](#compilation-and-loading)
  - [Prerequisites](#prerequisites)
  - [Linux / macOS](#linux--macos)
  - [Windows](#windows)
  - [Loading the Extension](#loading-the-extension)
- [Usage](#usage)
  - [Syntax](#syntax)
  - [Arguments](#arguments)
- [Examples](#examples)
  - [SQL Example](#sql-example)
  - [Python Example](#python-example)
- [Limitations](#limitations)

## How It Works

The code uses the Newton-Raphson method to find the discount rate at which the Net Present Value (NPV) of a series of cash flows equals zero. The function is implemented as an aggregate function and a window function, allowing it to be used in both contexts. The aggregate function computes the XIRR for all rows in a group, while the window function computes it for a specific subset of rows defined by a window frame.

### Key Features:

-   **Window Function:** Calculates XIRR for a moving window of data.
-   **Date Handling:** Automatically parses dates in common formats (`YYYY-MM-DD`, `DD/MM/YYYY`, `MM/DD/YYYY`, `YYYY.MM.DD`). Dates are converted to a numeric representation (days since January 1, 1900), with calculations using an average year length of `365.25` days to account for leap years.
-   **Optimization:** Uses the Newton-Raphson method with multiple starting points to increase the chances of finding a correct result.

## Compilation and Loading

To use this extension, you first need to compile it into a shared library.

### Prerequisites

Before compiling, you must have the SQLite development libraries installed, which provide the necessary header files (like `sqlite3.h`).

-   **On Debian/Ubuntu:** `sudo apt-get install sqlite3-dev`
-   **On Fedora/CentOS:** `sudo dnf install sqlite-devel`
-   **On macOS (with Homebrew):** `brew install sqlite`
-   **On Windows:** Ensure your compiler (like MinGW-w64) has access to the SQLite source code or pre-compiled headers.

### Linux / macOS

```bash
gcc -shared -fPIC -o sqlite-xirr-extension.so sqlite-xirr-extension.c -lm
```

### Windows

To compile on Windows, you can use a compiler like MinGW-w64 (GCC).

```bash
gcc -shared -o sqlite-xirr-extension.dll sqlite-xirr-extension.c -lm
```

### Loading the Extension

Once compiled, you can load the extension in your SQLite session:

```sql
-- On Linux/macOS
.load ./sqlite-xirr-extension.so

-- On Windows
.load ./sqlite-xirr-extension.dll
```

## Usage

The `xirr` function is available as an aggregate function and a window function. For convenience, the function is registered under two names: `xirr` (lowercase) and `XIRR` (uppercase).

### Syntax

The function accepts 2, 3, or 4 arguments.

```sql
-- Basic syntax (2 arguments)
XIRR(date, cash_flow)

-- With optional arguments
XIRR(date, cash_flow, [current_value], [initial_guesses])
```

As a window function:

```sql
XIRR(date, cash_flow, ...) OVER (PARTITION BY ... ORDER BY ...)
```

### Arguments

1.  `date` (text, **required**): The date of the cash flow.
2.  `cash_flow` (numeric, **required**): The value of the cash flow (positive for inflows, negative for outflows).
3.  `current_value` (numeric, *optional*): The total market value of the investment on the date of the transaction. This is treated as an additional cash flow at the end of the period, simulating a hypothetical sale of all assets. If omitted, it defaults to `0`.
4.  `initial_guesses` (text, *optional*): A string containing a custom list of starting rates for the Newton-Raphson algorithm, separated by a pipe (`|`). For example: `'0.1'` or `'0.1|0.25|-0.05'`. This allows you to guide the calculation if the default starting points fail to find a solution. By default, the function uses `'0.1|0.05|0.2|-0.1|0.01|0.5'` as the initial guesses list.

The function intelligently distinguishes between `current_value` and `initial_guesses` based on their data type. A numeric value will always be treated as `current_value`, and a text value will be treated as `initial_guesses`. You can provide one, both, or neither of these optional arguments.

## Examples

### SQL Example

Let's assume we have a `transactions` table with the following data:

| portfolio_id | transaction_date | amount | portfolio_value |
|--------------|------------------|--------|-----------------|
| 1            | 2023-01-01       | -1000  | 1000            |
| 1            | 2023-03-15       | -500   | 1600            |
| 1            | 2023-06-01       | 200    | 1900            |
| 1            | 2023-08-20       | 100    | 2100            |
| 2            | 2023-01-01       | -2000  | 2000            |
| 2            | 2023-04-01       | 500    | 2500            |
| 2            | 2023-07-01       | 150    | 3000            |

You can calculate the XIRR as follows:

#### Case 1: Basic Calculation (2 arguments)

This calculates the XIRR using only the transaction amounts. The final `current_value` is implicitly `0`.

```sql
SELECT
  portfolio_id,
  XIRR(transaction_date, amount) AS xirr_rate
FROM transactions
GROUP BY portfolio_id;
```

#### Case 2: With `current_value` (3 arguments)

This is the most common use case, where the running `portfolio_value` is included as the final cash flow in each calculation.

```sql
SELECT
  portfolio_id,
  XIRR(transaction_date, amount, portfolio_value) AS xirr_rate
FROM transactions
GROUP BY portfolio_id;
```

#### Case 3: With Custom Guesses (3 arguments)

If the calculation struggles to find a result, you can provide your own starting points.

```sql
SELECT
  portfolio_id,
  XIRR(transaction_date, amount, '0.1|0.05|-0.1') AS xirr_rate
FROM transactions
GROUP BY portfolio_id;
```

#### Case 4: With `current_value` and Custom Guesses (4 arguments)

You can provide both optional arguments. The function will identify them by type (numeric vs. text).

```sql
SELECT
  portfolio_id,
  XIRR(transaction_date, amount, portfolio_value, '0.15|0.2') AS xirr_rate
FROM transactions
GROUP BY portfolio_id;
```

#### Window Function Example

The function works seamlessly in a window context, calculating the XIRR for each row based on the expanding window of data. Arguments are the same as for the aggregate function.

```sql
SELECT
  portfolio_id,
  transaction_date,
  XIRR(transaction_date, amount, portfolio_value) OVER (
    PARTITION BY portfolio_id
    ORDER BY transaction_date
    ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
  ) AS xirr_rate
FROM transactions;
```

### Python Example

You can also load this extension directly in Python using the built-in `sqlite3` module.

```python
import sqlite3

# Path to your compiled extension
# Use './sqlite-xirr-extension.so' on Linux/macOS
# Use './sqlite-xirr-extension.dll' on Windows
EXTENSION_PATH = './sqlite-xirr-extension.so' 

# 1. Connect to your database
con = sqlite3.connect(':memory:') # Or a file path like 'my_database.db'

# 2. Enable extension loading (disabled by default for security)
con.enable_load_extension(True)

# 3. Load the extension
con.load_extension(EXTENSION_PATH)

print("Extension loaded successfully.")

# 4. You can now use the XIRR function in your queries
# (Example setup)
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

# Execute a query with the XIRR aggregate function (with current_value)
query = """
SELECT
  portfolio_id,
  XIRR(transaction_date, amount, portfolio_value) AS xirr_rate
FROM transactions
GROUP BY portfolio_id;
"""

print("\n--- XIRR with current_value ---")
for row in cur.execute(query):
    print(row)

# Execute a query with custom guesses
query = """
SELECT
  portfolio_id,
  XIRR(transaction_date, amount, portfolio_value, '0.1|0.2') AS xirr_rate
FROM transactions
GROUP BY portfolio_id;
"""

print("\n--- XIRR with custom guesses ---")
for row in cur.execute(query):
    print(row)

# 5. Close the connection
con.close()
```

## Limitations

-   **Date Range:** The function supports dates from **January 1, 1900** to **December 31, 9999**.
-   **Date Format:** Supported formats are `YYYY-MM-DD`, `DD/MM/YYYY`, `MM/DD/YYYY`, and `YYYY.MM.DD`. Other formats will not work.
-   **XIRR Result Range:** The calculated XIRR value is clamped to the range of **-0.9999** (-99.99%) to **10.0** (1000%). If the result falls outside this range, it may not be found.
-   **Prerequisite:** To calculate XIRR, at least one positive and one negative cash flow are required. Otherwise, the function will return `NULL`.
-   **Convergence:** In some cases, the Newton-Raphson algorithm may fail to converge to a solution. In such situations, the function will also return `NULL`.