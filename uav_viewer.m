classdef uav_viewer < handle
    %
    %    Create spacecraft animation
    %
    %--------------------------------
    properties
        body_handle
    	Vertices
    	Faces
    	facecolors
        plot_initialized
        t
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = uav_viewer
            self.body_handle = [];
            beta = [0, 0, 0]; % angle of attack
            [self.Vertices, self.Faces, self.facecolors] = self.define_spacecraft(beta);
            self.plot_initialized = 0;
            self.t = 0;
        end
        %---------------------------
        function self = update(self, time, state)
            if self.plot_initialized==0
                self = self.drawBody(state(1), state(2), state(3),...
                                   state(7), state(8), state(9),...
                                   state(13), state(14), state(15));
                xlabel('$\rho_e$ (m)', 'Interpreter', 'latex')
                ylabel('$\rho_n$ (m)', 'Interpreter', 'latex')
                zlabel('$h$ (m)', 'Interpreter', 'latex')
                axis([-1 + state(2), 1 + state(2),...
                      -1 + state(1), 1 + state(1),...
                      -1 - state(3), 1 - state(3)]);
                axis('equal')
                grid on
                self.plot_initialized = 1;
                title(gca, sprintf("Time: %.2f sec", self.t), "FontSize", 24);
            else
                self.t = time;
                self=self.drawBody(state(1), state(2), state(3),... 
                                   state(7), state(8), state(9),...
                                   state(13), state(14), state(15));
                axis([-1 + state(2), 1 + state(2),...
                      -1 + state(1), 1 + state(1),...
                      -1 - state(3), 1 - state(3)]);
                title(gca, sprintf("Time: %.2f sec", self.t), "FontSize", 24);
            end
        end
        %---------------------------
        function self = drawBody(self, pn, pe, pd, phi, theta, psi, alpha1, alpha2, alpha3)
            alpha = [alpha1, alpha2, alpha3];
            [self.Vertices, self.Faces, self.facecolors] = self.define_spacecraft(alpha);
            vertices = self.rotate(self.Vertices, phi, theta, psi);   % rotate rigid body  
            vertices = self.translate(vertices, pn, pe, pd);     % translate after rotation
            % transform vertices from NED to ENU (for matlab rendering)
            R = [...
                0, 1, 0;...
                1, 0, 0;...
                0, 0, -1;...
                ];
            vertices = R*vertices;
            if isempty(self.body_handle)
                self.body_handle = patch('Vertices', vertices', 'Faces', self.Faces,...
                                             'FaceVertexCData',self.facecolors,...
                                             'FaceColor','flat');
            else
                set(self.body_handle,'Vertices',vertices','Faces',self.Faces);
                drawnow
            end
        end 
        %---------------------------
        function pts = rotate(self, pts, phi, theta, psi)
            % define rotation matrix (right handed)
            R_roll = [...
                        1, 0, 0;...
                        0, cos(phi), sin(phi);...
                        0, -sin(phi), cos(phi)];
            R_pitch = [...
                        cos(theta), 0, -sin(theta);...
                        0, 1, 0;...
                        sin(theta), 0, cos(theta)];
            R_yaw = [...
                        cos(psi), sin(psi), 0;...
                        -sin(psi), cos(psi), 0;...
                        0, 0, 1];
            R = R_roll*R_pitch*R_yaw;   % inertial to body
            R = R';  % body to inertial
            % rotate vertices
            pts = R*pts;
        end
        %---------------------------
        % translate vertices by pn, pe, pd
        function pts = translate(self, pts, pn, pe, pd)
            pts = pts + repmat([pn;pe;pd],1,size(pts,2));
        end
        %---------------------------
        function [V, F, colors] = define_spacecraft(self, alpha)
            % Define the vertices (physical location of vertices in XYZ, inverted)
            V = [...
                  0.0000     0.0000    -0.0685;... % wing1 base top
                  0.0000     0.0000     0.0685;... % wing1 base bottom
                  0.4700     0.0000    -0.0515;... % wing1 end top
                  0.4700     0.0000     0.0585;... % wing1 end bottom
                  0.0000     0.0000     0.0000;... % wing2 base top
                  0.0000     0.0000     0.1370;... % wing2 base bottom
                 -0.2350     0.4070     0.0170;... % wing2 end top
                 -0.2350     0.4070     0.1270;... % wing2 end bottom
                  0.0000     0.0000     0.0000;... % wing3 base top
                  0.0000     0.0000     0.1370;... % wing3 base bottom
                 -0.2350    -0.4070     0.0170;... % wing3 end top
                 -0.2350    -0.4070     0.1270;... % wing3 end bottom
                  0.0750    -0.0750    -0.0750;... % box front top right
                  0.0750     0.0750    -0.0750;... % box front top left
                  0.0750    -0.0750     0.0750;... % box front bottom right
                  0.0750     0.0750     0.0750;... % box front bottom left
                 -0.0750    -0.0750    -0.0750;... % box back top right
                 -0.0750     0.0750    -0.0750;... % box back top left
                 -0.0750    -0.0750     0.0750;... % box back bottom right
                 -0.0750     0.0750     0.0750;... % box back bottom left
            ]';
            
            % Add Wing tilt
            V_link = V(:, 1:4);
            V(:, 1:4) = C(3, pi) * C(1, alpha(1)) * V_link;
            V(:, 5:8) = C(3, -pi/3) * C(1, alpha(2)) * V_link;
            V(:, 9:12) = C(3, pi/3) * C(1, alpha(3)) * V_link;

            % define faces as a list of vertices numbered above
            F = [...
                     1,  2,   4,   3;...  % wing1
                     5,  6,   8,   7;...  % wing2
                     9, 10,  12,  11;...  % wing3 
                    13, 14,  16,  15;...  % front
                    13, 15,  19,  17;...  % right
                    17, 18,  20,  19;...  % back
                    14, 18,  20,  16;...  % left
                    13, 14,  18,  17;...  % top
                    15, 16,  20,  19;...  % bottom
                    ];

            % define colors for each face    
            myred = [1, 0, 0]; 
            mygreen = [0, 1, 0];
            myblue = [0, 0, 1];
            myyellow = [1, 1, 0];
            mycyan = [0, 1, 1];
            myall = [1, 0, 1];

            colors = [...
                myred;... % wing1
                mygreen;...   % wing2
                myblue;...   % wing3
                myred;... % front
                myyellow;...   % right
                myyellow;...   % back
                myyellow;... % left
                mycyan;...   % top
                myall;...   % bottom                
                ];
        end
    end
end