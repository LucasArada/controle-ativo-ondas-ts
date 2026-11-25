function [Dx, Dy, F] = initDerivatives(X,Y)

method = 'SL6O3';
filter_alpha = 0.49;

[centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS] = finiteDifferenceCoefficients(method);
[filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS] = spatialFilterCoefficients(filter_alpha,true);

[LHSx, RHSx] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, {1}, {length(X)}, length(X));
[LHSy, RHSy] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, {1}, {length(Y)}, length(Y));

LHSx = applyMetric(LHSx,X);
LHSy = applyMetric(LHSy,Y);

LHSx = LHSx{1};
RHSx = RHSx{1};
LHSy = LHSy{1};
RHSy = RHSy{1};

[fLHSx, fRHSx] = makeMatricesEachDirection(filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS, {1}, {length(X)}, length(X));
[fLHSy, fRHSy] = makeMatricesEachDirection(filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS, {1}, {length(Y)}, length(Y));

fLHSx = fLHSx{1};
fRHSx = fRHSx{1};
fLHSy = fLHSy{1};
fRHSy = fRHSy{1};

Dx = @(U)(LHSx\(RHSx*U));
Dy = @(U)((LHSy\(RHSy*U'))');

F = @(U)((fLHSy\(fRHSy*(fLHSx\(fRHSx*U))'))');

end

function [LHS, RHS] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStarts, derivEnds, n)

    nTypes = length(derivStarts);
    
    LHS = cell(1,nTypes);
    RHS = cell(1,nTypes);
    
    if n == 1 % If single point in this direction, set derivative to zero and filter to one
        if centeredStencilRHS(1) == 0
            for i = 1:nTypes
                LHS{i} = 1;
                RHS{i} = 0;
            end
        else
            for i = 1:nTypes
                LHS{i} = 1;
                RHS{i} = 1;
            end
        end
        return
    end
    
    % Create the base stencils
    LHS_base = diag(centeredStencilLHS(1)*ones(n,1));
    RHS_base = diag(centeredStencilRHS(1)*ones(n,1));
    
    if centeredStencilRHS(1) == 0 % If the RHS centered stencil has center 0, its left side should have the signal inverted as this is the stencil for derivatives, otherwise, its a stencil for filters
        invertStencil = -1;
    else
        invertStencil = 1;
    end
        
    for i = 1:length(centeredStencilLHS)-1
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1+i,n));
        LHS_base(inds) = centeredStencilLHS(i+1);
        
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1-i,n));
        LHS_base(inds) = centeredStencilLHS(i+1);
    end
    
    for i = 1:length(centeredStencilRHS)-1
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1+i,n));
        RHS_base(inds) = centeredStencilRHS(i+1);
        
        inds = sub2ind([n n], 1:n, 1+mod((1:n)-1-i,n));
        RHS_base(inds) = invertStencil*centeredStencilRHS(i+1);
    end
    
    % Add startings and endings
    
    [mLHS, nLHS] = size(decenteredStencilLHS);
    [mRHS, nRHS] = size(decenteredStencilRHS);
    
    for i = 1:nTypes
        LHS_temp = LHS_base;
        RHS_temp = RHS_base;
        
        for j = 1:length(derivStarts{i})
            ind_start = derivStarts{i}(j);
            
            LHS_temp(ind_start:ind_start+mLHS-1,:) = 0;
            RHS_temp(ind_start:ind_start+mRHS-1,:) = 0;
            
            LHS_temp(ind_start:ind_start+mLHS-1,ind_start:ind_start+nLHS-1) = decenteredStencilLHS;
            RHS_temp(ind_start:ind_start+mRHS-1,ind_start:ind_start+nRHS-1) = decenteredStencilRHS;
            
        end
        
        for j = 1:length(derivEnds{i})
            ind_end = derivEnds{i}(j);
            
            LHS_temp(ind_end-mLHS+1:ind_end,:) = 0;
            RHS_temp(ind_end-mRHS+1:ind_end,:) = 0;
            
            LHS_temp(ind_end-mLHS+1:ind_end,ind_end-nLHS+1:ind_end) = decenteredStencilLHS(end:-1:1,end:-1:1);
            RHS_temp(ind_end-mRHS+1:ind_end,ind_end-nRHS+1:ind_end) = invertStencil*decenteredStencilRHS(end:-1:1,end:-1:1);
            
        end
        
        LHS{i} = sparse(LHS_temp);
        RHS{i} = sparse(RHS_temp);
        
    end
end

function LHS = applyMetric(LHS,X)

    n = length(X);
    if n == 1
        return
    end
    [centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS] = finiteDifferenceCoefficients('SL4');
    [LHS_temp, RHS_temp] = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, {1}, {n}, n);
    
    dXdEta = LHS_temp{1}\(RHS_temp{1}*X');
    
    nTypes = length(LHS);
    for i = 1:nTypes
        LHS{i} = LHS{i}.*dXdEta';
    end

end

function [centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS] = finiteDifferenceCoefficients(method)

% This functions contains the coefficients for each finite differences method. New methods can be defiend here if needed.

% Centered Stencils start at the center point and move to the right hand side
% Uncentered stencils are as they should be in the left hand side top of the matrix

switch method
        
	case 'SL6'

		w = 1.8;

		% Center
		matriz_a = [1	   1	      1			 -2;
					1	   4	      9			 -6;
					1	   16	      81		 -10;
					sin(w) sin(2*w)/2 sin(3*w)/3 -2*w*cos(w)];

		matriz_b = [1; 0; 0; w];

		coeffs = matriz_a\matriz_b;

		a     = coeffs(1)/2;
        b     = coeffs(2)/4;
		c	  = coeffs(3)/6;	        
		alpha = coeffs(4);

		centeredStencilLHS = [1 alpha];
        centeredStencilRHS = [0 a b c];
		
		%
		%
		% Border
		% First stage
		matriz_a = [1	1	1	1	1		1		0;
					0	1	2	3	4		5		-1;
					0	1	4	9	16		25		-2;
					0	1	8	27	64		125		-3;
					0	1	16	81	256		625		-4;
					0	1	32	243	1024	3125	-5;
					0	1	64	729	4096	15625	-6];
		
		matriz_b = [0; 1; 0; 0; 0; 0; 0]; %[0; -1; 0; 0; 0; 0; 0];

		coeffs = matriz_a\matriz_b;

		a     = coeffs(1);
        b     = coeffs(2);
		c	  = coeffs(3);
		d     = coeffs(4);
        e     = coeffs(5);
		f	  = coeffs(6);	        
		alpha = coeffs(7);
		
		matriz_lhs_aux = [1 alpha 0 0];
		matriz_rhs_aux = [a b c d e f];

		% Second stage
		matriz_a = [-1	2;
					-1	6];

		matriz_b = [-1; 0];
		
		coeffs = matriz_a\matriz_b;

		a     = coeffs(1)/2;
		alpha = coeffs(2);

		matriz_lhs_aux2 = [alpha 1 alpha 0];
		matriz_rhs_aux2 = [-a 0 a 0 0 0];

		% Third stage (same as centered SL4)

		
		matriz_a = [1      0          -2/3       ;
                    0      1          -4/3       ;
                    sin(w) sin(2*w)/2 -2*w*cos(w)];
        
        matriz_b = [4/3; -1/3; w];
        
        coeffs = matriz_a\matriz_b;
        
        a     = coeffs(1)/2;
        b     = coeffs(2)/4;
        alpha = coeffs(3);
        
        matriz_lhs_aux3 = [0 alpha 1 alpha];
        matriz_rhs_aux3 = [-b -a 0 a b 0];

		
		decenteredStencilLHS = [matriz_lhs_aux;
								matriz_lhs_aux2;
								matriz_lhs_aux3];         
        
        decenteredStencilRHS = [matriz_rhs_aux;
								matriz_rhs_aux2;
								matriz_rhs_aux3];

	case 'SL6O3'
		% Inner-Scheme: SL6 - Optimized, alpha determined minimizing the integral of dispersion
		%                     from w = [0,wf], i.e., alpha = f(integral(Edisp)|_0^wf)), wf = pi/2;
		% P3: SL4 - Optimized following the procedure above, wf = pi/2;
		% P2: C4 - Pade Scheme
		% P1: C3 -  Optimized, max. dissipation error Ediss(w=pi) = pi/2; (Edisp(w=pi/2)~6%)
		coeffs =[0.392465753424658, 1.565410958904110, 0.237260273972603, -0.017739726027397];
		alpha = coeffs(1); a = coeffs(2); b = coeffs(3); c = coeffs(4);
		centeredStencilLHS  = [1 , alpha];
		centeredStencilRHS  = [0 , a/2, b/4, c/6];

		coeffs_P3 = [0.350978473581213, 1.567318982387476, 0.134637964774951];
		alpha_P3 = coeffs_P3(1); a_P3 = coeffs_P3(2); b_P3 = coeffs_P3(3);

		decenteredStencilLHS = [    1  , (3*pi + 40)/(3*pi + 8),   0   ,     0     ;...
								  1/4 ,            1          ,  1/4   ,     0     ;...
								   0  ,        alpha_P3       ,   1   ,  alpha_P3];
		decenteredStencilRHS = [ -(13*pi + 56)/(2*(3*pi + 8)), (15*pi + 8)/(2*(3*pi + 8)), -(3*pi - 56)/(2*(3*pi + 8)), (pi - 8)/(2*(3*pi + 8)), 0;...
									 -3/4  ,    0    ,   3/4  ,     0         0   ;...
								  -b_P3/4  , -a_P3/2 ,    0   ,  a_P3/2 ,  b_P3/4 ];
										  
    case 'SL4' % 4th order spectral-like compact derivatives
        
        w = 1.8;
		           
        matriz_a = [1      0          -2/3       ;
                    0      1          -4/3       ;
                    sin(w) sin(2*w)/2 -2*w*cos(w)];
        
        matriz_b = [4/3; -1/3; w];
        
        coeffs = matriz_a\matriz_b;
        
        a     = coeffs(1)/2;
        b     = coeffs(2)/4;
        alpha = coeffs(3);
        
        centeredStencilLHS = [1 alpha];
        centeredStencilRHS = [0 a b];
        
        
        decenteredStencilLHS = [1   3 0  ;
                                1/4 1 1/4];
        
        decenteredStencilRHS = [-17/6 3/2 3/2 -1/6;
                                -3/4  0   3/4 0   ];
            
        
    case 'EX2' % 2nd order explicit finite differences
    
        centeredStencilLHS = 1;
        centeredStencilRHS = [0 1/2];
        
        decenteredStencilLHS = 1;
        decenteredStencilRHS = [-3/2 2 -1/2];
        
    case 'EX4' % 4nd order explicit finite differences
    
        centeredStencilLHS = 1;
        centeredStencilRHS = [0 2/3 -1/12];
        
        decenteredStencilLHS = 1;
        decenteredStencilRHS = [-25/12 4 -3 4/3 -1/4
                                -1/2   0 1/2 0    0 ];
        
    otherwise
        error('Finite differences method not implemented: %s. Check finiteDifferenceCoefficients file for available methods',method)
end

end

function [filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS] = spatialFilterCoefficients(alpha,filterBorders)

% This function outputs the coefficients for the 10th order spatial filter from Gaitonde 1998 for a given alpha

filterStencilLHS = [1 alpha];

filterStencilRHS = [(193+126*alpha)/256 (105+302*alpha)/512 (-15+30*alpha)/128 (45-90*alpha)/1024 (-5+10*alpha)/512 (1-2*alpha)/1024];

if islogical(filterBorders)
	if filterBorders
		filterBorders = 'decentered';
	else
		filterBorders = 'off';
	end
end

switch filterBorders
	case 'decentered'
		filterDecenteredStencilLHS = [1 alpha];
		
		a_bound( 1,5) = (-1+2*alpha)/1024;
		a_bound( 2,5) = (5-10*alpha)/512;
		a_bound( 3,5) = (-45+90*alpha)/1024;
		a_bound( 4,5) = (15+98*alpha)/128;
		a_bound( 5,5) = (407+210*alpha)/512;
		a_bound( 6,5) = (63+130*alpha)/256;
		a_bound( 7,5) = (-105+210*alpha)/512;
		a_bound( 8,5) = (15-30*alpha)/128;
		a_bound( 9,5) = (-45+90*alpha)/1024;
		a_bound(10,5) = (5-10*alpha)/512;
		a_bound(11,5) = (-1+2*alpha)/1024;

		a_bound( 1,4) = (1-2*alpha)/1024;
		a_bound( 2,4) = (-5+10*alpha)/512;
		a_bound( 3,4) = (45+934*alpha)/1024;
		a_bound( 4,4) = (113+30*alpha)/128;
		a_bound( 5,4) = (105+302*alpha)/512;
		a_bound( 6,4) = (-63+126*alpha)/256;
		a_bound( 7,4) = (105-210*alpha)/512;
		a_bound( 8,4) = (-15+30*alpha)/128;
		a_bound( 9,4) = (45-90*alpha)/1024;
		a_bound(10,4) = (-5+10*alpha)/512;
		a_bound(11,4) = (1-2*alpha)/1024;

		a_bound( 1,3) = (-1+2*alpha)/1024;
		a_bound( 2,3) = (5+502*alpha)/512;
		a_bound( 3,3) = (979+90*alpha)/1024;
		a_bound( 4,3) = (15+98*alpha)/128;
		a_bound( 5,3) = (-105+210*alpha)/512;
		a_bound( 6,3) = (63-126*alpha)/256;
		a_bound( 7,3) = (-105+210*alpha)/512;
		a_bound( 8,3) = (15-30*alpha)/128;
		a_bound( 9,3) = (-45+90*alpha)/1024;
		a_bound(10,3) = (5-10*alpha)/512;
		a_bound(11,3) = (-1+2*alpha)/1024;

		a_bound( 1,2) = (1+1022*alpha)/1024;
		a_bound( 2,2) = (507+10*alpha)/512;
		a_bound( 3,2) = (45+934*alpha)/1024;
		a_bound( 4,2) = (-15+30*alpha)/128;
		a_bound( 5,2) = (105-210*alpha)/512;
		a_bound( 6,2) = (-63+126*alpha)/256;
		a_bound( 7,2) = (105-210*alpha)/512;
		a_bound( 8,2) = (-15+30*alpha)/128;
		a_bound( 9,2) = (45-90*alpha)/1024;
		a_bound(10,2) = (-5+10*alpha)/512;
		a_bound(11,2) = (1-2*alpha)/1024;

		a_bound( 1,1) = (1023+1*alpha)/1024;
		a_bound( 2,1) = (5+507*alpha)/512;
		a_bound( 3,1) = (-45+45*alpha)/1024;
		a_bound( 4,1) = (15-15*alpha)/128;
		a_bound( 5,1) = (-105+105*alpha)/512;
		a_bound( 6,1) = (63-63*alpha)/256;
		a_bound( 7,1) = (-105+105*alpha)/512;
		a_bound( 8,1) = (15-15*alpha)/128;
		a_bound( 9,1) = (-45+45*alpha)/1024;
		a_bound(10,1) = (5-5*alpha)/512;
		a_bound(11,1) = (-1+1*alpha)/1024;
		
		filterDecenteredStencilRHS = a_bound';
		
	case 'reducedOrder'
        filterDecenteredStencilLHS = 1;
        
		F2 = [1/2+alpha, 1/2+alpha]./[1 2];
		F4 = [5/8+3/4*alpha, 1/2+alpha, -1/8+1/4*alpha]./[1 2 2];
		F6 = [11/16+5/8*alpha, 15/32+17/16*alpha, -3/16+3/8*alpha, 1/32-1/16*alpha]./[1 2 2 2];
		F8 = [93/128+70/128*alpha, 7/16+18/16*alpha, -7/32+14/32*alpha, 1/16-1/8*alpha, -1/128+1/64*alpha]./[1 2 2 2 2];
		
		filterDecenteredStencilRHS = zeros(5,9);
		
		filterDecenteredStencilRHS(1,1) = 1;
		filterDecenteredStencilRHS(2,1:3) = F2([2 1 2]);
		filterDecenteredStencilRHS(3,1:5) = F4([3 2 1 2 3]);
		filterDecenteredStencilRHS(4,1:7) = F6([4 3 2 1 2 3 4]);
		filterDecenteredStencilRHS(5,1:9) = F8([5 4 3 2 1 2 3 4 5]);
		
	case 'off'
		filterDecenteredStencilLHS = diag(ones(1,5));
		filterDecenteredStencilRHS = diag(ones(1,5));
end

end

