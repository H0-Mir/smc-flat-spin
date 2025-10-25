# Aircraft Flat-Spin Recovery Using Sliding Mode Control

This repository provides a MATLAB implementation of the control algorithm proposed in the paper  
**“Aircraft Flat-Spin Recovery Using Sliding-Mode Based Attitude and Altitude Control.”**  
The simulation demonstrates the recovery of an aircraft from a flat-spin condition using a nonlinear sliding-mode–based control strategy.

---

## 🧩 Control Law

The control input is defined as:

$$
\tilde{u} = -
\begin{bmatrix}  
	g_{\delta_e}^{y_i} & g_{\delta_a}^{y_i} & g_{\delta_r}^{y_i} & g_{\delta_t}^{y_i}  
\end{bmatrix}^{-1}
\left[
f^{y_i} + k_i(\dot{y}_i - \dot{y}_{d_i}) - \ddot{y}_{d_i} + \lambda_i \left| \tilde{s}_i \right|^{a_i} \text{sgn}(\tilde{s}_i)
\right]
$$

---

## 📘 Stability Analysis

The Lyapunov candidate function is defined as:

$$
\tilde{V}_{i} = \frac{1}{2} \tilde{s}_{i}^\top \tilde{s}_{i}
\implies
\dot{\tilde{V}}_{i} = \tilde{s}_{i}^\top(-\lambda_{i}|\tilde{s}_{i}|^{a_{i}} \text{sgn}(\tilde{s}_{i}))
= -\lambda_{i}|\tilde{s}_{i}|^{a_{i}+1}
$$

---

## 🚀 Running the Simulation

1. Clone or download this repository.  
2. Open MATLAB and navigate to the project directory.  
3. Run the following command in the MATLAB Command Window:

```matlab
main
```
---
## 📚 References

> Salahudden, Dwivedi, V. S., Dwivedi, P. N., Giri, D. K., & Ghosh, A. K. (2021).  
> *Aircraft Flat-Spin Recovery Using Sliding-Mode Based Attitude and Altitude Control.*  
> *Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering*, **235**(8), 924–936.  
> [SAGE Publications, London, UK]  
> [https://doi.org/10.1177/0954410020964674](https://doi.org/10.1177/0954410020964674)

---



