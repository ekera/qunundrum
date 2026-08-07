# The <code>solve_linear_distribution_shor</code> executable

## Synopsis
```console
Synopsis: mpirun solve_linear_distribution_shor
   [ -t-bound <t-bound> ] [ -cofactor-bound <cofactor-bound> ]
      [ -instances <instances> ]
         <distribution> { <distribution> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs for the order $r$ using Shor's original post-processing algorithm based on continued fraction expansion.

In total $10^3$ problem instances are considered by default to gather statistics; this may be changed via the <code>-instances</code> flag.

The results are written to the console and to a log file <code>logs/solve-linear-shor-YYYYMMDD-HHmmss±ZZZZ-XXXXXXXX.txt</code> where <code>YYYYMMDD-HHmmss</code> is the date and time when the executable was started, <code>±ZZZZ</code> is the timezone offset from UTC in hours and minutes (HHmm), and <code>XXXXXXXX</code> is a random 32-bit integer in hexadecimal. This prevents concurrently running executables from writing to the same log file.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

### Optional command line arguments
Flag specifying the search bound on $t$ (defaults to $2^8$):
- <code>-t-bound \<t-bound\></code> sets the search bound on $|t|$ to <code>\<t-bound\></code>; all integers $j + t$ are solved for the order $r$

Flag specifying the search bound on the cofactor (defaults to $2^{16}$):
- <code>-cofactor-bound \<cofactor-bound\></code> sets the bound on the cofactor between the order $r$ and the denominator to <code>\<cofactor-bound\></code>

Flag specifying the number of problem instances to solve (defaults to $10^3$):
- <code>-instances \<instances\></code> sets the number of problem instances to <code>\<instances\></code>

## Interpreting the output
The log file is on the format
```
# Log file: solve-linear-shor-20240229-014632+0100-6c1b4d72.txt

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-1.txt
# Bounds: (t: 256, cofactor: 65536)
# Timestamp: 2024-02-29 01:46:32 CET
m: 2048 s: 1 n: 1 -- success: 999 -- fail: 1 (0) -- prepare:    20.329 ms solve:    20.699 ms [    9.937,  5111.765]

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-1.txt
# Bounds: (t: 0, cofactor: 1)
# Timestamp: 2024-02-29 01:48:36 CET
m: 2048 s: 1 n: 1 -- success: 220 -- fail: 780 (0) -- prepare:     0.949 ms solve:    11.039 ms [    9.886,    12.118]
```
where we find $m$, $s$ or $\ell$, $n$ — #success — #fail — prep-time — solve-time, and where
- $m$ is the bit length of the order $r$,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated; otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- #success is the number of problem instances that were successfully solved,
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors,
- prep-time is the average time in ms required to setup the problem instances,
- solve-time is the average [min, max] time in ms required to solve the problem instances.
