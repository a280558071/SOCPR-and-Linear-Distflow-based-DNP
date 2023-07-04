# SOCPR-and-Linear-Distflow-based-DNP
Two examples for distribution network planning (DNP) method based on Second-Order cone programming (SOCP) relaxation and Linear Distflow are included here. 

The Objective is to minimize the investment cost on distribution lines and the only operation cost, value of lost load(VOLL). 

The modeling is completed by YALMIP, which is a toolbox in MATLAB to do optimization modeling*, and solved by GUROBI. 
Therefore, you should have MATLAB, YALMIP and GUROBI installed in your PC to make these codes work.

"LinearDistFlow_DNP.m" is DNP based on Linear Distflow Optimal Power Flow (OPF) Model in [1] and [2], but also described in a compact form in Section III in [3].

"SOCPROPF_DNP.m" is DNP based on SOCP Relaxation OPF Model in [3].

* *See https://yalmip.github.io/ for more information.
* [1] M. E. Baran and F. F. Wu, “Optimal capacitor placement on radial distribution systems,” IEEE Trans. Power Delivery, vol. 4, no. 1, pp. 725–734,1989.
* [2] M. E. Baran and F. F. Wu, “Optimal sizing of capacitors placed on a radial distribution system,” IEEE Trans. Power Delivery, vol. 4, no. 1, pp. 735–743, 1989
* [3] L Gan, N Li, U Topcu, SH Low, Exact convex relaxation of optimal power flow in radial networks, IEEE Transactions on Automatic Control 60 (1), 72-87, 2014
* Update in 2022/07/26 is important, previous version "SOCPROPF_DNP.m" could not produce optimal results due to a mistake in modeling.
