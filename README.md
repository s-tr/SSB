# SSB - a library for computing Singapore Savings Bonds rates

## Usage

```python
>>> import SSB

>>> # Average rates for 1, 2, 5, 10 years
>>> benchmark = [0.038, 0.035, 0.032, 0.033])
>>> rates = SSB.ssb(benchmark)
>>> print([round(100*rate, 2) for rate in rates])
[3.18, 3.18, 3.18, 3.18, 3.18, 3.18, 3.3, 3.5, 3.62, 3.62]
```


