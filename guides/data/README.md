# Data for the Guides

[`macro3_Q.csv`](macro3_Q.csv) holds three quarterly US series from 1960Q1 to 2019Q4, 240
quarters, with each quarter dated by its first day:

| Column | Series | Source (FRED) |
|---|---|---|
| `GDPC1` | real GDP growth, annualized, in percent: 400 times the log difference of real GDP | GDPC1 |
| `PCECTPI` | PCE inflation, annualized, in percent: 400 times the log difference of the PCE chain-type price index | PCECTPI |
| `UNRATE` | the unemployment rate, in percent, averaged over the three months of the quarter | UNRATE |

The values are the first 240 rows of the GDP growth, PCE inflation and unemployment columns of
`chapter12/macro4_Q.csv` in the code repository of Chan (2027), which dates a quarter by its
last month. They agree with the FRED series as published on 15 October 2026 to within 5e-9 for
GDP growth and PCE inflation, and to the four decimals of the file for unemployment. The sample
ends before the pandemic, as in the forecasting application of that book.

The guides of the series Forecasting with a VAR fit a VAR to these data and simulate from it.

[`soe_Q.csv`](soe_Q.csv) holds four US and four Australian quarterly series from 1984Q1 to
2019Q4, 144 quarters, dated in the same way:

| Column | Series | Source |
|---|---|---|
| `GDPC1` | US real GDP growth, annualized, in percent: 400 times the log difference of real GDP | FRED, GDPC1 |
| `CPIAUCSL` | US CPI inflation, annualized, in percent: 400 times the log difference of the CPI averaged over the three months of the quarter | FRED, CPIAUCSL |
| `FEDFUNDS` | the federal funds rate, in percent, averaged over the three months of the quarter | FRED, FEDFUNDS |
| `EBP` | the US excess bond premium, in percentage points, averaged over the three months of the quarter | Federal Reserve Board, [`ebp_csv.csv`](https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv) |
| `NGDPRSAXDCAUQ` | Australian real GDP growth, annualized, in percent | FRED, NGDPRSAXDCAUQ |
| `AUSCPIALLQINMEI` | Australian CPI inflation, annualized, in percent | FRED, AUSCPIALLQINMEI |
| `IR3TIB01AUQ156N` | the Australian 3-month interbank rate, in percent | FRED, IR3TIB01AUQ156N |
| `CCUSMA02AUQ618N` | the depreciation of the Australian dollar against the US dollar, annualized, in percent: 400 times the log difference of the Australian dollars per US dollar | FRED, CCUSMA02AUQ618N |

The FRED series were retrieved on 15 November 2026 and the excess bond premium on 16 November
2026. The premium is the part of the corporate bond spread of Gilchrist and Zakrajšek (2012) that
expected defaults do not explain. The Federal Reserve Board publishes the monthly series with the
note of Favara, Gilchrist, Lewis and Zakrajšek (2016). The sample starts after the Australian
dollar was floated in December 1983 and ends before the pandemic. Two quarters of Australian CPI
inflation stand out: Medicare lowered the CPI in 1984Q1, and the goods and services tax raised it
in 2000Q3. The Australian CPI is not seasonally adjusted.

The guides of the series VARs for a Small Open Economy fit a VAR to these data and simulate from
it.

## References

Chan, J. C. C. (2027). *Bayesian Macroeconometrics: Methods and Applications*. Chapman &
Hall/CRC. [Code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics)

Favara, G., Gilchrist, S., Lewis, K. F. and Zakrajšek, E. (2016). Updating the Recession Risk and
the Excess Bond Premium. *FEDS Notes*, Board of Governors of the Federal Reserve System, October
2016. [doi:10.17016/2380-7172.1836](https://doi.org/10.17016/2380-7172.1836)

Gilchrist, S. and Zakrajšek, E. (2012). Credit Spreads and Business Cycle Fluctuations. *American
Economic Review*, 102(4): 1692-1720.
[doi:10.1257/aer.102.4.1692](https://doi.org/10.1257/aer.102.4.1692)

U.S. Bureau of Economic Analysis, Real Gross Domestic Product [GDPC1] and Personal Consumption
Expenditures: Chain-type Price Index [PCECTPI]; U.S. Bureau of Labor Statistics, Unemployment Rate
[UNRATE]; retrieved from FRED, Federal Reserve Bank of St. Louis.

U.S. Bureau of Economic Analysis, Real Gross Domestic Product [GDPC1]; U.S. Bureau of Labor
Statistics, Consumer Price Index for All Urban Consumers: All Items in U.S. City Average
[CPIAUCSL]; Board of Governors of the Federal Reserve System (US), Federal Funds Effective Rate
[FEDFUNDS]; International Monetary Fund, Real Gross Domestic Product (GDP) for Australia
[NGDPRSAXDCAUQ]; Organization for Economic Co-operation and Development, Consumer Price Indices
(CPIs, HICPs), COICOP 1999: Consumer Price Index: Total for Australia [AUSCPIALLQINMEI], Interest
Rates: 3-Month or 90-Day Rates and Yields: Interbank Rates: Total for Australia
[IR3TIB01AUQ156N], and Currency Conversions: US Dollar Exchange Rate: Average of Daily Rates:
National Currency: USD for Australia [CCUSMA02AUQ618N]; retrieved from FRED, Federal Reserve Bank
of St. Louis.
