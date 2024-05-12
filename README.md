# SSB - a library for computing Singapore Savings Bonds rates

## Usage

```python
>>> import SSB
>>> # Average rates for 1, 2, 5, 10 years
>>> benchmark = [3.80, 3.50, 3.20, 3.30]
>>> rates = SSB.ssb(benchmark)
>>> print(rates)
[3.18, 3.18, 3.18, 3.18, 3.18, 3.18, 3.3, 3.5, 3.62, 3.62]
```


