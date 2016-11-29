The 'A' and 'B' results are for a proton-electron plasma with beta_||p = 1.0 and 0.01 respectively.
For the other plasma parameters we have equal proton and electron temperatures, with isotropy perpendicular and parallel to the magnetic field.

For completeness, I have attached the solutions for the Alfven, Fast, Entropy, and Slow modes, respectively modes 1 through 4.
The wavevector scan is split across two files, 
*k_1_1_1000_100.mode* 
scans from (kperp,kpar) rho_p = (1.E-3,1.E-3) to (1.,1.E-1)
while 
*k_1000_100_10000_215.mode* 
scans from (kperp,kpar) rho_p = (1.,1.E-1) to (1.E1,2.15)

The data is formatted as

kperp,kpar,beta,vtp/c
!1, 2, 3, 4
omega,gamma
!5, 6
Bx  By   Bz         Ex,     Ey     Ez
!7,8 9,10 11,12 : 13,14 15,16 17,18 : 
Upx,    Upy    Upz     Uex    Uey   Uez
!19,20 21,22 23,24 : 25,26 27,28 29,30 :
np       ne
!31,32 33,34
proton parameters
           !35-40
electron parameters
           !41-46

MVA of a single mode will be the square of 7-10 divided by the square of 11 and 12.
