### Build time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| python | mathoptmodel-batch | 0.2 ms | 1.6 ms | 36.7 ms |
| python | mathoptmodel-row | 17.3 ms | 165.5 ms | 2.12 s |
| python | pyomo | 95.4 ms | 1.30 s | 14.71 s |
| python | pyoptinterface-row | 126.8 ms | 1.27 s | 13.26 s |
| python | pyoptinterface-matrix | 145.9 ms | 1.49 s | 15.98 s |
| java | mathoptmodel-batch | 0.6 ms | 2.1 ms | 37.9 ms |
| java | mathoptmodel-row | 3.8 ms | 18.4 ms | 246.8 ms |
| java | ortools-modelbuilder | 9.1 ms | 107.0 ms | 1.01 s |
| java | ortools-mpsolver | 7.6 ms | 118.4 ms | 1.90 s |

### LP write time

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| python | mathoptmodel-batch | 18.0 ms | 168.6 ms | 2.47 s |
| python | mathoptmodel-row | 18.2 ms | 171.7 ms | 2.20 s |
| python | pyomo | 91.2 ms | 1.39 s | 18.87 s |
| python | pyoptinterface-row | 18.7 ms | 224.5 ms | 3.81 s |
| python | pyoptinterface-matrix | 18.4 ms | 216.7 ms | 3.73 s |
| java | mathoptmodel-batch | 13.2 ms | 150.7 ms | 2.06 s |
| java | mathoptmodel-row | 13.5 ms | 151.0 ms | 2.04 s |
| java | ortools-modelbuilder | 9.0 ms | 104.6 ms | 1.76 s |
| java | ortools-mpsolver | 16.0 ms | 280.6 ms | 4.40 s |

### Peak memory

| Language | Tool | n = 10,000 | n = 100,000 | n = 1,000,000 |
|---|---|---|---|---|
| python | mathoptmodel-batch | 30 MB | 75 MB | 422 MB |
| python | mathoptmodel-row | 45 MB | 112 MB | 835 MB |
| python | pyomo | 70 MB | 283 MB | 2115 MB |
| python | pyoptinterface-row | 61 MB | 153 MB | 1185 MB |
| python | pyoptinterface-matrix | 63 MB | 126 MB | 722 MB |
| java | mathoptmodel-batch | 80 MB | 147 MB | 767 MB |
| java | mathoptmodel-row | 92 MB | 184 MB | 944 MB |
| java | ortools-modelbuilder | 110 MB | 250 MB | 1356 MB |
| java | ortools-mpsolver | 115 MB | 333 MB | 2236 MB |

### Models written (read back with mathoptmodel)

| Language | Tool | Check |
|---|---|---|
| python | mathoptmodel-batch | ok |
| python | mathoptmodel-row | ok |
| python | pyomo | ok |
| python | pyoptinterface-row | ok |
| python | pyoptinterface-matrix | ok |
| java | mathoptmodel-batch | ok |
| java | mathoptmodel-row | ok |
| java | ortools-modelbuilder | ok |
| java | ortools-mpsolver | ok |
