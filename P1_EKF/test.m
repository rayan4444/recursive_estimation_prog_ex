function partial_diff

syms sx 
syms sy
syms rho
syms phi
syms u 
syms C_da
syms C_dh
syms C_w


%f= cos(phi)*(tanh(u)-C_dh*(sx^2+sy^2))- C_da*(sx - C_w*cos(rho))*sqrt((sx-C_w*cos(rho))^2+(sy-C_w*sin(rho))^2);
f= sin(phi)*(tanh(u)-C_dh*(sx^2+sy^2))- C_da*(sy - C_w*sin(rho))*sqrt((sx-C_w*cos(rho))^2+(sy-C_w*sin(rho))^2);
diff(f, rho)
end
