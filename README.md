# iterated-Arnoldi-Tikhonov
Code for the iterated Arnoldi-Tikhonov method presented in [1].

1. Select a Problem for Testing: Choose from the problems studied in [1] - phillips [2], baart [3], or blur [3].
2. Set Noise Level: Specify the relative value of noise ($\xi$ in [1]).
3. Determine Arnoldi Steps: Set the number of Arnoldi steps ($\ell$ in [1]). The decomposition is computed using the "Arnoldi" function as described in [1].
4. Compute Equation RHS: The program calculates the right-hand side of equation (14) from [1].
5. Iterated Tikhonov Steps: Decide on the number of steps for the iterated Tikhonov process. Using $i=1$ implements the Arnoldi-Tikhonov method presented in [2].
6. Condition Checking: The program verifies the condition of equation (15) in [1]. (Tip: To apply the condition for the parameter choice of the Arnoldi-Tikhonov method proposed in [2], add the commented line `E=3*E`).
7. Parameter Estimation: The "parameter" function estimates the method's parameter ($\alpha$ in [1]) using bisection.
8. Solution and Error Calculation: The program computes the solution and its relative error.
9. Visualization Option: The user is prompted to choose whether to plot the computed solution alongside the exact solution in one dimension.

[1] D. Bianchi, M. Donatelli, D. Furchì, L. Reichel. “Convergence analysis and parameter estimation for the iterated Arnoldi-Tikhonov method”, arXiv identifier: 2311.11823

[2] R. Ramlau, L. Reichel. “Error estimates for Arnoldi-Tikhonov”, In: Inverse Problems 35 (2019), p. 055002.

[3] P. C. Hansen. “Regularization tools version 4.0 for Matlab 7.3”, In: Numerical Algorithms 46 (2007), pp. 189-194.
