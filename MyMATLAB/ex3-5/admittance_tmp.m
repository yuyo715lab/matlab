%define impedance
r11 = -30i;
r12 = 0.08+0.4i;
r13 = 0.12+0.5i;
r22 = -34i;
r21 = r12;
r23 = 0.1+0.4i;
r33 = -29i;
r31 = r13;
r32 = r23;
r34 = 0.3i;

%winding ratio
n=1.1;

%direct admittance
Y11 = 1/r12 + 1/r11 + 1/r13;
Y22 = 1/r22 + 1/r21 + 1/r23;
Y33_= 1/r33 + 1/r31 + 1/r32 + 1/r34;
Y33 = Y33_ + (n^2-1)*1/r34;
Y44 = 1/r34;

%transfer admittance
Y12 = -1/r12;
Y13 = -1/r13;
Y14 = 0;
Y21 = Y12;
Y23 = -1/r23;
Y24 = 0;
Y31 = Y13;
Y32 = Y23;
Y34_=-1/r34;
Y34 = -(n-1)*1/r34 + Y34_;
Y41 = Y14;
Y42 = Y24;
Y43 = Y34;


%admittance matrix
Y = [Y11 Y12 Y13 Y14
  Y21 Y22 Y23 Y24
  Y31 Y32 Y33 Y34
  Y41 Y42 Y43 Y44]