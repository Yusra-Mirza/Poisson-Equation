/******************************************************
* PROJECT: Numerical Solution of 2D Poisson Equation
*
* DESCRIPTION:
* This program solves the 2D Poisson equation using
* three iterative numerical methods:
*
* 1. Jacobi Method
* 2. Gauss-Seidel Method
* 3. Successive Over-Relaxation (SOR)
*
* Users can choose the method at runtime.
* The potential values are saved into a CSV file for
* visualization or post-processing.
*
* EQUATION FORM:
* ∇²φ = -ρ
*
* APPLICATIONS:
* - Electrostatics
* - Heat distribution
* - Fluid flow
* - Computational physics
*
* AUTHOR: Arpit, Shivani, Aryan, Sachin, Shubh, Mustufa
******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NX 50
#define NY 50
#define MAX_ITER 10000
#define TOL 1e-6
#define OMEGA 1.5 // Relaxation factor for SOR
// ---------------- Jacobi Method ----------------
void solveJacobi(double phi[NX][NY], double rho[NX][NY], double h) {
double phi_new[NX][NY] = {0};
double error;
for (int iter = 0; iter < MAX_ITER; iter++) {
error = 0.0;
// Update using only values from OLD iteration
for (int i = 1; i < NX - 1; i++) {
for (int j = 1; j < NY - 1; j++) {
phi_new[i][j] = 0.25 * (
phi[i+1][j] + phi[i-1][j] +
phi[i][j+1] + phi[i][j-1] -
h * h * rho[i][j]
);
double diff = fabs(phi_new[i][j] - phi[i][j]);
if (diff > error) error = diff;
}
}
// Copy new values back
for (int i = 1; i < NX - 1; i++)
for (int j = 1; j < NY - 1; j++)
phi[i][j] = phi_new[i][j];
// Convergence check
if (error < TOL) {
printf("Jacobi converged after %d iterations, error = %e\n", iter, error);
return;
}
}
}
// ---------------- Gauss-Seidel Method ----------------
void solveGaussSeidel(double phi[NX][NY], double rho[NX][NY], double h) {
double error;
for (int iter = 0; iter < MAX_ITER; iter++) {
error = 0.0;
// Updates use latest values (in-place)
for (int i = 1; i < NX - 1; i++) {
for (int j = 1; j < NY - 1; j++) {
double old = phi[i][j];
phi[i][j] = 0.25 * (
phi[i+1][j] + phi[i-1][j] +
phi[i][j+1] + phi[i][j-1] -
h * h * rho[i][j]
);
double diff = fabs(phi[i][j] - old);
if (diff > error) error = diff;
}
}
if (error < TOL) {
printf("Gauss-Seidel converged after %d iterations, error = %e\n", iter, error);
return;
}
}
}
// ---------------- SOR Method ----------------
void solveSOR(double phi[NX][NY], double rho[NX][NY], double h) {
double error;
for (int iter = 0; iter < MAX_ITER; iter++) {
error = 0.0;
for (int i = 1; i < NX - 1; i++) {
for (int j = 1; j < NY - 1; j++) {
double old = phi[i][j];
double gs_val = 0.25 * (
phi[i+1][j] + phi[i-1][j] +
phi[i][j+1] + phi[i][j-1] -
h * h * rho[i][j]
);
// SOR update with relaxation
phi[i][j] = (1 - OMEGA) * old + OMEGA * gs_val;
double diff = fabs(phi[i][j] - old);
if (diff > error) error = diff;
}
}
if (error < TOL) {
printf("SOR converged after %d iterations, error = %e\n", iter, error);
return;
}
}
}
// ---------------- Main Program ----------------
int main() {
double phi[NX][NY] = {0}; // Potential grid
double rho[NX][NY] = {0}; // Source field
double h = 1.0 / (NX - 1); // Grid spacing
// Place a point source at grid center
int cx = NX / 2, cy = NY / 2;
rho[cx][cy] = 1.0;
// Get user choice
int choice;
printf("Select method:\n");
printf("1. Jacobi Method\n");
printf("2. Gauss-Seidel Method\n");
printf("3. SOR Method (with Relaxation)\n");
printf("Enter choice: ");
scanf("%d", &choice);
// Run selected method
switch (choice) {
case 1: solveJacobi(phi, rho, h); break;
case 2: solveGaussSeidel(phi, rho, h); break;
case 3: solveSOR(phi, rho, h); break;
default:
printf("Invalid choice!\n");
return 0;
}
// Save numerical solution to CSV
FILE *fp = fopen("phi_output.csv", "w");
for (int i = 0; i < NX; i++) {
for (int j = 0; j < NY; j++) {
fprintf(fp, "%f", phi[i][j]);
if (j < NY - 1) fprintf(fp, ",");
}
fprintf(fp, "\n");
}
fclose(fp);
printf("Data saved to phi_output.csv\n");
return 0;
}