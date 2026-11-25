### Build time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| python | mathoptmodel-batch | 0.2 ms | 3.2 ms | 41.1 ms |
| python | mathoptmodel-row | 17.0 ms | 172.1 ms | 1.73 s |
| python | pyomo | 94.8 ms | 1.31 s | 14.72 s |
| python | pyomo-linear-expression | 58.8 ms | 860.6 ms | 9.94 s |
| python | pyoptinterface-row | 127.6 ms | 1.29 s | 13.50 s |
| python | pyoptinterface-matrix | 145.7 ms | 1.50 s | 16.22 s |
| java | mathoptmodel-batch | 0.5 ms | 4.2 ms | 50.5 ms |
| java | mathoptmodel-row | 6.6 ms | 18.7 ms | 287.4 ms |
| java | ortools-modelbuilder | 10.0 ms | 87.4 ms | 1.00 s |
| java | ortools-modelbuilder-helper | 4.8 ms | 36.7 ms | 262.0 ms |
| java | ortools-mpsolver | 8.6 ms | 123.1 ms | 1.94 s |
| java | ortools-mpsolver-proto | 10.0 ms | 113.6 ms | 1.52 s |

### LP write time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| python | mathoptmodel-batch | 1.9 ms | 16.7 ms | 182.4 ms |
| python | mathoptmodel-row | 2.2 ms | 16.2 ms | 173.2 ms |
| python | pyomo | 93.4 ms | 1.43 s | 19.07 s |
| python | pyomo-linear-expression | 90.8 ms | 1.39 s | 19.32 s |
| python | pyoptinterface-row | 18.7 ms | 242.0 ms | 4.06 s |
| python | pyoptinterface-matrix | 18.5 ms | 223.2 ms | 3.91 s |
| java | mathoptmodel-batch | 1.8 ms | 18.8 ms | 186.0 ms |
| java | mathoptmodel-row | 1.9 ms | 16.5 ms | 177.6 ms |
| java | ortools-modelbuilder | 10.1 ms | 106.5 ms | 1.83 s |
| java | ortools-modelbuilder-helper | 8.3 ms | 104.4 ms | 1.74 s |
| java | ortools-mpsolver | 15.6 ms | 278.7 ms | 4.24 s |
| java | ortools-mpsolver-proto | 16.3 ms | 254.1 ms | 4.24 s |

### Peak memory

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| python | mathoptmodel-batch | 29 MB | 46 MB | 213 MB |
| python | mathoptmodel-row | 45 MB | 85 MB | 600 MB |
| python | pyomo | 71 MB | 283 MB | 2115 MB |
| python | pyomo-linear-expression | 71 MB | 280 MB | 2153 MB |
| python | pyoptinterface-row | 71 MB | 153 MB | 1184 MB |
| python | pyoptinterface-matrix | 71 MB | 127 MB | 722 MB |
| java | mathoptmodel-batch | 77 MB | 103 MB | 268 MB |
| java | mathoptmodel-row | 91 MB | 149 MB | 420 MB |
| java | ortools-modelbuilder | 110 MB | 239 MB | 1280 MB |
| java | ortools-modelbuilder-helper | 95 MB | 172 MB | 876 MB |
| java | ortools-mpsolver | 115 MB | 332 MB | 2193 MB |
| java | ortools-mpsolver-proto | 135 MB | 380 MB | 2339 MB |

### Models written (read back with mathoptmodel)

| Language | Tool | Check |
|---|---|---|
| python | mathoptmodel-batch | ok |
| python | mathoptmodel-row | ok |
| python | pyomo | ok |
| python | pyomo-linear-expression | ok |
| python | pyoptinterface-row | ok |
| python | pyoptinterface-matrix | ok |
| java | mathoptmodel-batch | ok |
| java | mathoptmodel-row | ok |
| java | ortools-modelbuilder | ok |
| java | ortools-modelbuilder-helper | ok |
| java | ortools-mpsolver | ok |
| java | ortools-mpsolver-proto | ok |
