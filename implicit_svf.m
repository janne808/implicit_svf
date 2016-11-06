%  (C) 2016 Janne Heikkarainen <janne.heikkarainen@tut.fi>
%
%  All rights reserved.
%
%  This file is part of implicit_svf.m, a trapezoidal implicit integrator
%  for a state-variable filter
%
%  implicit_svf.m is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  implicit_svf.m is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with implicit_svf.m. If not, see <http://www.gnu.org/licenses/>.
%

clear all;
close all;

% samplerate
fs=44100;

% input signal
input=wavread('saw.wav');
input=0.6*input;
tt=1:length(input)-1;

% output signal
output=zeros(1,length(tt));

% cutoff modulation
ff=linspace(1.5,0.0,length(tt));

% filter state
hp=0;
bp=0;
lp=0;

% cutoff variable
f=0.25;

% q variable
q=0.1;

for tt=tt
    % cutoff frequency modulation
    f=ff(tt);
    
    % integrate filter state
    A_n=bp+0.5*f*tanh(-lp-q*bp+input(tt));
    aa=lp;
    bb=bp;
    cc=input(tt+1);
    xx=bp;
    bp_1=bp;
    % iterate to find root with newton's method
    for mm=1:32
        dd=cc-aa-q*xx-0.5*f*(tanh(bb)+tanh(xx));
        xx2=xx-(xx-0.5*f*tanh(dd)-A_n)/...
                 (1-0.5*f*(q-0.5*f*(tanh(xx)^2-1))*(tanh(dd)^2-1));
        %disp(['step: ' num2str(mm) ' xx: ' num2str(xx)...
        %      ' dx:' num2str(xx2-xx)]);
        if(abs(xx2-xx)<1e-15)
            xx=xx2;
            break;
        end
        xx=xx2;
    end
    bp=xx;
    lp=lp+0.5*f*(tanh(bp_1)+tanh(bp));
    
    % output
    output(tt)=lp;
end

% plot
plot(output);

% wav write
wavwrite(0.25*output,fs,'output.wav');

