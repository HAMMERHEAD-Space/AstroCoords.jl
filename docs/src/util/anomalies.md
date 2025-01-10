# Anomalies

There is support in the package to switch between the various anomalies -- True, Eccentric, and Mean along with a Kepler Solver based on a Newton iterative method. These may be moved out to a separate package later to better support and experiment with various Kepler Solver methodologies. 

![image](../assets/anomalies.png)
*Comparison of Anomalies [1]*

# Usage

To convert between the anomalies use a simple function call works.

```julia
M = 2π * rand()
e = rand()

E = meanAnomaly2EccentricAnomaly(M, e)
```

## References
[1]: https://www.bogan.ca/orbits/kepler/e_anomly.html
