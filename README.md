Wavelet_Multiscale_Correlation
==============================

The Wavelet Multi-scale Correlation (WMC) uses wavelet coefficients to calculatecorrelations between time series at different time scales. It has been used in widerange of fields like finance, economics, and climate studies [1, 3]. 

Fourier transform does not represent abrupt changes efficiently because it is not alocal function in time. Wavelets are well localized in time and frequency can be usedfor simultaneous time and frequency analysis. Using a Discrete Wavelet Transform(DWT), a time series can be decomposed into time and frequency components. The decomposed components are called wavelet coefficients. Here we focus on themethodology of using these wavelet co-efficient to obtain correlations at differenttime scales. For more detailed introduction of wavelets and its applications, I can recommend reading the books [3, 2]. 

A type of DWT called Maximal Overlap Discrete Wavelet Transform (MODWT) gives same number of wavelet co-efficient at all the time scales and hence is used in correlation analysis. Consider two discrete time series $`X`$ and $`Y`$ of length $`T`$. Using MODWT, each of the time series can be decomposed into $`J = 1,2,...\log_{2}T`$ time scales. Each of these time scales will have $`T`$ discrete wavelet coefficients, each represented as $`d^{j,t}_X`$ and $`d^{j,t}_Y`$, where j is the time scale and $`t=1,2,...,T`$ is the discretetime step. Variance and co-variance of these time series can be then calculated as:

src="https://latex.codecogs.com/svg.latex?Var_{X}^{j}\equiv\frac{1}{T_j}\sum_{t=M_j-1}^{T-1}\left[d^{j,t}_{X}\right]^2" title="Var_{X}^{j}\equiv\frac{1}{T_j}\sum_{t=M_j-1}^{T-1}\left[d^{j,t}_{X}\right]^2" /></a>

where $Var_{X}^{j}$ and $Var_{Y}^{j}$ are the variances at $j^{th}$ scale of the time series $X$ and $Y$ respectively. $COV_{XY}^{j}$ is the co-variance between $X$ and $Y$ at $j^{th}$ scale. \mbox{$T_j=T-M_j-1$} stands for the number of wavelet coefficients unaffected by the boundary, with $M_j = (2^j-1)(M-1)$ and $M$ is the length of the wavelet filter used. The \acf{wmc} at $j^{th}$ scale can then be found using the equation:
\begin{equation}
    \rho_{XY}^{j}\equiv \frac{COV_{XY}^{j}}{\sqrt{Var_{Y}^{j}Var_{X}^{j}}},
\end{equation}
The values of $\rho_{XY}^{j}$ can be between -1 and 1, with the extreme values indicating complete correlation or inverse correlation, and 0 indicating no correlation.
