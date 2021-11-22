#=
Solves DSGE using second order approximation perturbation methods

Mauricio M Tejada
Department of Economics, Alberto Hurtado University
March, 2020
---

Based on:
[1] S. Schmitt-Grohe, M. Uribe, Solving dynamic general equilibrium models using a second-order approximation to the policy function, Journal of Economic Dynamics and Control 28 (2004) 755–775.
[2] Caraiani, Petre,  Introduction to Quantitative Macroeconomics Using Julia, 2019, Academic Press  Elsevier

---

Functions:
J, H  = compute_gradient_hessian(model_eqs, xs)
Gs,Hs,Gx,Hx = solve_first_order_approx(xs,J,nx,ny,ne)
Gs,Hs,Gx,Hx,Gxx,Hxx,Gss,Hss = solve_second_order_approx(xs,J,H,ETA,nx,ny,ne)

=#
using ForwardDiff
using LinearAlgebra

function compute_gradient_hessian(model_eqs, xs)

    neqs = length(model_eqs);
    nvars = length(xs)
    J = Array{Float64}(undef, neqs, nvars)
    H = Array{Float64}(undef, nvars*nvars, neqs)

    for i in 1:neqs
        Jeq = ForwardDiff.gradient(model_eqs[i], xs);
        J[i,:] = Jeq'

        Heq = ForwardDiff.hessian(model_eqs[i], xs);
        H[:,i] = vec(Heq)
    end

    return J, H

end

function solve_first_order_approx(xs,J,nx,ny,ne)

    Hs=xs[1:nx]; 
    Gs=xs[nx+1:nx+ny];

    # Order 1 
    J0=J[:,1:nx+ny]; 
    J1=-J[:,nx+ny+1:end];

    # Complex Schur Decomposition
    qz = schur(J0,J1); 

    # Pick non-explosive (stable) eigenvalues
    slt = abs.(diag(qz.T)).<abs.(diag(qz.S))
    nk=sum(slt);

    # Reorder the system with stable eigs in upper-left
    qzo = ordschur(qz, slt)

    # Split up the results appropriately
    z21 = qzo.Z[nk+1:end,1:nk];
    z11 = qzo.Z[1:nk,1:nk];

    s11 = qzo.S[1:nk,1:nk];
    t11 = qzo.T[1:nk,1:nk];

    if rank(z11)<nx; 
        error("Invertibility condition violated") 
    end

   #compute Gx and Hx
    z11i = z11\Matrix{Float64}(I, nk, nk);
    Gx = z21*z11i;  
    Hx = z11*(s11\t11)*z11i;

    tol=1e-6; 
    if maximum(maximum(imag(Gx)))<tol 
        Gx=real(Gx); 
    end 
    if maximum(maximum(imag(Gx)))<tol 
        Hx=real(Hx); 
    end

    return Gs,Hs,Gx,Hx

end

function solve_second_order_approx(xs,J,H,ETA,nx,ny,ne)

    Hs=xs[1:nx]; 
    Gs=xs[nx+1:nx+ny];

    # Order 1 
    J0=J[:,1:nx+ny]; 
    J1=-J[:,nx+ny+1:end];

    # Complex Schur Decomposition
    qz = schur(J0,J1); 

    # Pick non-explosive (stable) eigenvalues
    slt = abs.(diag(qz.T)).<abs.(diag(qz.S))
    nk=sum(slt);

    # Reorder the system with stable eigs in upper-left
    qzo = ordschur(qz, slt)

    # Split up the results appropriately
    z21 = qzo.Z[nk+1:end,1:nk];
    z11 = qzo.Z[1:nk,1:nk];

    s11 = qzo.S[1:nk,1:nk];
    t11 = qzo.T[1:nk,1:nk];

    if rank(z11)<nx; 
        error("Invertibility condition violated") 
    end

   #compute Gx and Hx
    z11i = z11\Matrix{Float64}(I, nk, nk);
    Gx = z21*z11i;  
    Hx = z11*(s11\t11)*z11i;

    tol=1e-6; 
    if maximum(maximum(imag(Gx)))<tol 
        Gx=real(Gx); 
    end 
    if maximum(maximum(imag(Gx)))<tol 
        Hx=real(Hx); 
    end

    # Computes Gxx and Hxx
    Zx = [Hx;Gx * Hx;Matrix{Float64}(I, nx, nx);Gx];
    Jxp= J[:,1:nx]; Jyp= J[:,nx+1:nx+ny];
    Jx = J[:,nx+ny+1:2 * nx+ny];
    Jy = J[:,2 * nx+ny+1:2 * (nx+ny)]; 
    XX1 = [kron((Jxp+Jyp * Gx),Matrix{Float64}(I,nx * nx,nx * nx)) kron(Jyp,kron(Hx',Hx'))+kron(Jy,Matrix{Float64}(I,nx * nx,nx * nx))]; 
    XX0 = -kron(Zx',Zx') * H; 
    XX0 = vec(XX0); 
    HGXX= \(XX1,XX0); 
    Hxx = HGXX[1:nx * nx * nx]; 
    
    if maximum(maximum(imag(Hxx)))<tol 
        Hxx=real(Hxx); 
    end 
    
    Gxx = HGXX[nx * nx * nx+1:end]; 
    
    if maximum(maximum(imag(Gxx)))<tol 
        Gxx=real(Gxx); 
    end

    # Computes Gss and Hss 
    SS0= 0.0; 
    
    for i=1:ne; 
        Zs =[ETA[:,i];
        Gx * ETA[:,i];zeros(nx,1);zeros(ny,1)]; 
        TEMP0 = -kron(Zs',Zs') * H; 
        TEMP0 = vec(TEMP0); 
        TEMP1 = -Jyp * kron(Matrix{Float64}(I,ny,ny),kron(ETA[:,i]',ETA[:,i]')) * Gxx; 
        TEMP1 = vec(TEMP1); 
        SS0 = SS0.+TEMP0.+TEMP1; 
    end

    SS1 = [Jxp+Jyp * Gx (Jyp+Jy)]; 
    HGSS = \(SS1,SS0); 
    Hss = HGSS[1:nx];

    Gss = HGSS[nx+1:nx+ny]; 
    
    if maximum(maximum(imag(Hss)))<tol 
        Hss=real(Hss); 
    end 
    
    if maximum(maximum(imag(Gss)))<tol 
        Gss=real(Gss); 
    end

    Hxx = reshape(Hxx,nx * nx,nx)'; 
    Gxx = reshape(Gxx,nx * nx,ny)';

    return Gs,Hs,Gx,Hx,Gxx,Hxx,Gss,Hss

end

function fir_first_order(Hx,Gx,num_periodos,shock,η,σ)

    # Horizonte
    horizonte = 50
    t = 1:horizonte

    # El shock inicial es de una desviación estándar.
    shock_inicial = η[:,shock].*σ   

    # Situación inicial: Estado estacionario
    ny, = size(Gx)
    nx, = size(Hx)
    xx = zeros(horizonte,nx);
    yy = zeros(horizonte,ny);

    for i = 1:horizonte
        xx[i,:] = (Hx^i)*shock_inicial
        yy[i,:] = Gx*(Hx^i)*shock_inicial
    end

    return xx, yy

end