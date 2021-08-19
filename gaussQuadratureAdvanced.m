% ............................................................. 

    function [weights,locations]=gaussQuadratureAdvanced(option)
    % Gauss quadrature for Q4 elements
    % option 'complete' (2x2)
    % option 'reduced'  (1x1)
    % locations: Gauss point locations
    % weights: Gauss point weights
        
    switch option
        case 'complete'
    
        locations=...
          [ -0.577350269189626 -0.577350269189626;
             0.577350269189626 -0.577350269189626;
             0.577350269189626  0.577350269189626;
            -0.577350269189626  0.577350269189626];
        weights=[ 1;1;1;1]; 
    
        case 'reduced'
        
        locations=[0 0];
        weights=[4];
        
        case 'advanced'
        locations=...
            [-0.774596669241483 -0.774596669241483;
             -0.774596669241483 0.0;
             -0.774596669241483 0.774596669241483;
             0. -0.774596669241483;
            0. 0.0;
            0. 0.774596669241483;
            0.774596669241483 -0.774596669241483;
            0.774596669241483 0.0;
            0.774596669241483 0.774596669241483];
        weights=...
            [0.555555555555556*0.555555555555556;
0.555555555555556*0.888888888888889;
0.555555555555556*0.555555555555556;
0.555555555555556*0.888888888888889;
0.888888888888889*0.888888888888889;
0.555555555555556*0.555555555555556;
0.555555555555556*0.555555555555556;
0.555555555555556*0.888888888888889;
0.555555555555556*0.555555555555556];
                
    end

    end  % end function gaussQuadrature
    