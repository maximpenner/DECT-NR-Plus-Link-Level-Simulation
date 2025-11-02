function [W, codebook_index_max] = Beamforming_W(N_TS, N_TX, codebook_index)

    q = 1i;

    % 6.3.4 and Table 6.3.4-1
    if N_TS == 1 && N_TX == 1
        
        W = 1;

        codebook_index_max = 0;
    
    elseif N_TS == 1 && N_TX == 2
        
        switch codebook_index
            case 0
                W = [1;0];
            case 1
                W = [0;1];
            case 2
                W = 1/sqrt(2)*[1;1];
            case 3
                W = 1/sqrt(2)*[1;-1];
            case 4
                W = 1/sqrt(2)*[1;q];
            case 5
                W = 1/sqrt(2)*[1;-q];
            otherwise
                error('Unknwon codebook index %d.', codebook_index);
        end

        codebook_index_max = 5;
        
    % Table 6.3.4-2
    elseif N_TS == 1 && N_TX == 4
        
        switch codebook_index
            case 0
                W = [1;0;0;0];
            case 1
                W = [0;1;0;0];
            case 2
                W = [0;0;1;0];
            case 3
                W = [0;0;0;1];
                
                case 4
                    W = 1/sqrt(2)*[1;0;1;0];
                case 5
                    W = 1/sqrt(2)*[1;0;-1;0];
                case 6
                    W = 1/sqrt(2)*[1;0;q;0];
                case 7
                    W = 1/sqrt(2)*[1;0;-q;0];
                case 8
                    W = 1/sqrt(2)*[0;1;0;1];
                case 9
                    W = 1/sqrt(2)*[0;1;0;-1];
                case 10
                    W = 1/sqrt(2)*[0;1;0;q];
                case 11
                    W = 1/sqrt(2)*[0;1;0;-q];
                
                    case 12
                        W = 1/2*[1;1;1;1];
                    case 13
                        W = 1/2*[1;1;q;q];
                    case 14
                        W = 1/2*[1;1;-1;-1];
                    case 15
                        W = 1/2*[1;1;-q;-q];
                    case 16
                        W = 1/2*[1;q;1;q];
                    case 17
                        W = 1/2*[1;q;q;-1];
                    case 18
                        W = 1/2*[1;q;-1;-q];
                    case 19
                        W = 1/2*[1;q;-q;1];
                    case 20
                        W = 1/2*[1;-1;1;-1];
                    case 21
                        W = 1/2*[1;-1;q;-q];
                    case 22
                        W = 1/2*[1;-1;-1;1];
                    case 23
                        W = 1/2*[1;-1;-q;q];
                    case 24
                        W = 1/2*[1;-q;1;-q];
                    case 25
                        W = 1/2*[1;-q;q;1];
                    case 26
                        W = 1/2*[1;-q;-1;q];
                    case 27
                        W = 1/2*[1;-q;-q;-1];
                
            otherwise
                error('Unknwon codebook index %d.', codebook_index);
        end

        codebook_index_max = 27;
        
    % Table 6.3.4-3
    elseif N_TS == 2 && N_TX == 2
        
        switch codebook_index
            case 0
                W = 1/sqrt(2)*[1 0; 0 1];
            case 1
                W = 1/2*[1 1; 1 -1];
            case 2
                W = 1/2*[1 1; q -q];                
                
            otherwise
                error('Unknown codebook index %d.', codebook_index);
        end

        codebook_index_max = 2;
        
    % Table 6.3.4-4
    elseif N_TS == 2 && N_TX == 4
        
        switch codebook_index
            case 0
                W = 1/sqrt(2)*[ 1  0;...
                                0  1;...
                                0  0;...
                                0  0];
            case 1
                W = 1/sqrt(2)*[ 1  0;...
                                0  0;...
                                0  1;...
                                0  0];
            case 2
                W = 1/sqrt(2)*[ 1  0;...
                                0  0;...
                                0  0;...
                                0  1];
            case 3
                W = 1/sqrt(2)*[ 0  0;...
                                1  0;...
                                0  1;...
                                0  0];
            case 4
                W = 1/sqrt(2)*[ 0  0;...
                                1  0;...
                                0  0;...
                                0  1];
            case 5
                W = 1/sqrt(2)*[ 0  0;...
                                0  0;...
                                1  0;...
                                0  1];
                            
                case 6
                    W = 1/2*[   1  0;...
                                0  1;...
                                1  0;...
                                0 -q];
                case 7
                    W = 1/2*[   1  0;...
                                0  1;...
                                1  0;...
                                0  q];
                case 8
                    W = 1/2*[   1  0;...
                                0  1;...
                               -q  0;...
                                0  1];
                case 9
                    W = 1/2*[   1  0;...
                                0  1;...
                               -q  0;...
                                0 -1];
                case 10
                    W = 1/2*[   1  0;...
                                0  1;...
                               -1  0;...
                                0 -q];
                case 11
                    W = 1/2*[   1  0;...
                                0  1;...
                               -1  0;...
                                0  q];
                case 12
                    W = 1/2*[   1  0;...
                                0  1;...
                                q  0;...
                                0  1];
                case 13
                    W = 1/2*[   1  0;...
                                0  1;...
                                q  0;...
                                0 -1];
                    case 14
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                            1  1;...
                                            1 -1;...
                                            1 -1];
                    case 15
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                            1  1;...
                                            q -q;...
                                            q -q];
                    case 16
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                            q  q;...
                                            1 -1;...
                                            q -q];
                    case 17
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                            q  q;...
                                            q -q;...
                                           -1  1];
                    case 18
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                           -1 -1;...
                                            1 -1;...
                                           -1  1];
                    case 19
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                           -1 -1;...
                                            q -q;...
                                           -q  q];
                    case 20
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                           -q -q;...
                                            1 -1;...
                                           -q  q];
                    case 21
                        W = 1/(2*sqrt(2))*[ 1  1;...
                                           -q -q;...
                                            q -q;...
                                            1 -1];
                
            otherwise
                error('Unknwon codebook index %d.', codebook_index);
        end

        codebook_index_max = 21;
        
    % Table 6.3.4-2
    elseif N_TS == 4 && N_TX == 4
        
        switch codebook_index
            case 0
                W = 1/2*[   1  0  0  0;...
                            0  1  0  0;...
                            0  0  1  0;...
                            0  0  0  1];
                case 1
                    W = 1/(2*sqrt(2))*[ 1  1  0  0;...
                                        0  0  1  1;...
                                        1 -1  0  0;...
                                        0  0  1 -1];
                case 2
                    W = 1/(2*sqrt(2))*[ 1  1  0  0;...
                                        0  0  1  1;...
                                        q -q  0  0;...
                                        0  0  q -q];
                    case 3
                        W = 1/4*[   1  1  1  1;...
                                    1 -1  1 -1;...
                                    1  1 -1 -1;...
                                    1 -1 -1  1];
                    case 4
                        W = 1/4*[   1  1  1  1;...
                                    1 -1  1 -1;...
                                    q  q -q -q;...
                                    q -q -q  q];
                
            otherwise
                error('Unknwon codebook index %d.', codebook_index);
        end

        codebook_index_max = 4;
        
    % Table 6.3.4-2
    elseif N_TS == 8 && N_TX == 8
        
        switch codebook_index
            case 0
                %W = 1/4*eye(8);             % error in TS
                W = 1/(2*sqrt(2))*eye(8);
                
            otherwise
                error('Unknwon codebook index %d.', codebook_index);
        end

        codebook_index_max = 1;
        
    else
        error('Beamforming configuration is unknown, N_TS = %d and N_TX = %d.', N_TS, N_TX);
    end
end
