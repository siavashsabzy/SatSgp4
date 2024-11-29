classdef SatSgp4

    properties
        satrec
        name
        iar80Ar
        rar80Ar
    end


    methods

        function obj = SatSgp4(longstr1, longstr2, nameString)
            nut80 = obj.getNut80();
            obj.iar80Ar = nut80(:,1:5);
            obj.rar80Ar = nut80(:,6:9);
            for i=1:106
                for j=1:4
                    obj.rar80Ar(i,j)= obj.rar80Ar(i,j) * 0.0001 * pi / (180*3600.0);
                end
            end
            obj.name = nameString;
            obj.satrec = obj.parseTle(longstr1, longstr2);
        end

        function eph = propagate(obj, julianDay)
            % output in ECI
            daysToGo = julianDay - obj.satrec.jDayEpoch;
            if daysToGo < 0

            end
            [~, r, v] = obj.propagateSingleTeme( obj.satrec, daysToGo*1440);
            [reci, veci] = obj.temeToEci( r, v, julianDay);
            eph = [reci', veci'] .*1000.0;
        end



    end
    methods (Access= private)
        function [satrec, r, v] = propagateSingleTeme(obj, satrec, tsince)

            % /* ------------------ set mathematical constants --------------- */
            twopi = 2.0 * pi;
            x2o3  = 2.0 / 3.0;
            % sgp4fix divisor for divide by zero check on inclination
            % the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
            % 1.5 e-12, so the threshold was changed to 1.5e-12 for consistancy
            temp4    =   1.5e-12;

            %  // sgp4fix identify constants and allow alternate values
            % no longer needed, pass in through satrec
            %global tumin mu radiusearthkm xke j2 j3 j4 j3oj2
            vkmpersec     = satrec.radiusearthkm * satrec.xke/60.0;

            % /* --------------------- clear sgp4 error flag ----------------- */
            satrec.t     = tsince;
            satrec.error = 0;
            mrt = 0.0;

            % /* ------- update for secular gravity and atmospheric drag ----- */
            xmdf    = satrec.mo + satrec.mdot * satrec.t;
            argpdf  = satrec.argpo + satrec.argpdot * satrec.t;
            nodedf  = satrec.nodeo + satrec.nodedot * satrec.t;
            argpm   = argpdf;
            mm      = xmdf;
            t2      = satrec.t * satrec.t;
            nodem   = nodedf + satrec.nodecf * t2;
            tempa   = 1.0 - satrec.cc1 * satrec.t;
            tempe   = satrec.bstar * satrec.cc4 * satrec.t;
            templ   = satrec.t2cof * t2;

            if (satrec.isimp ~= 1)
                delomg = satrec.omgcof * satrec.t;
                delm   = satrec.xmcof *...
                    ((1.0 + satrec.eta * cos(xmdf))^3 -...
                    satrec.delmo);
                temp   = delomg + delm;
                mm     = xmdf + temp;
                argpm  = argpdf - temp;
                t3     = t2 * satrec.t;
                t4     = t3 * satrec.t;
                tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 -...
                    satrec.d4 * t4;
                tempe  = tempe + satrec.bstar * satrec.cc5 * (sin(mm) -...
                    satrec.sinmao);
                templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +...
                    satrec.t * satrec.t5cof);
            end

            nm    = satrec.no;
            em    = satrec.ecco;
            inclm = satrec.inclo;
            if (satrec.method == 'd')
                tc = satrec.t;
                [satrec.atime,em,argpm,inclm,satrec.xli,mm,...
                    satrec.xni,nodem,dndt,nm] = obj.dspace(...
                    satrec.d2201,satrec.d2211,satrec.d3210,...
                    satrec.d3222,satrec.d4410,satrec.d4422,...
                    satrec.d5220,satrec.d5232,satrec.d5421,...
                    satrec.d5433,satrec.dedt,satrec.del1,...
                    satrec.del2,satrec.del3,satrec.didt,...
                    satrec.dmdt,satrec.dnodt,satrec.domdt,...
                    satrec.irez,satrec.argpo,satrec.argpdot,satrec.t,...
                    tc,satrec.gsto,satrec.xfact,satrec.xlamo,satrec.no,...
                    satrec.atime,em,argpm,inclm,satrec.xli,mm,...
                    satrec.xni,nodem,nm);
            end % // if method = d

            if (nm <= 0.0)
                %       fprintf(1,'# error nm %f\n', nm);
                satrec.error = 2;
            end
            am = (satrec.xke / nm)^x2o3 * tempa * tempa;
            nm = satrec.xke / am^1.5;
            em = em - tempe;

            % // fix tolerance for error recognition
            if ((em >= 1.0) || (em < -0.001) || (am < 0.95))
                %       fprintf(1,'# error em %f\n', em);
                satrec.error = 1;
            end
            %   sgp4fix change test condition for eccentricity
            if (em < 1.0e-6)
                em  = 1.0e-6;
            end
            mm     = mm + satrec.no * templ;
            xlm    = mm + argpm + nodem;
            emsq   = em * em;
            temp   = 1.0 - emsq;
            nodem  = rem(nodem, twopi);
            argpm  = rem(argpm, twopi);
            xlm    = rem(xlm, twopi);
            mm     = rem(xlm - argpm - nodem, twopi);

            % /* ----------------- compute extra mean quantities ------------- */
            sinim = sin(inclm);
            cosim = cos(inclm);

            % /* -------------------- add lunar-solar periodics -------------- */
            ep     = em;
            xincp  = inclm;
            argpp  = argpm;
            nodep  = nodem;
            mp     = mm;
            sinip  = sinim;
            cosip  = cosim;
            if (satrec.method == 'd')
                [ep,xincp,nodep,argpp,mp] = obj.dpper(...
                    satrec.e3,satrec.ee2,satrec.peo,...
                    satrec.pgho,satrec.pho,satrec.pinco,...
                    satrec.plo,satrec.se2,satrec.se3,...
                    satrec.sgh2,satrec.sgh3,satrec.sgh4,...
                    satrec.sh2,satrec.sh3,satrec.si2,...
                    satrec.si3,satrec.sl2,satrec.sl3,...
                    satrec.sl4,satrec.t,satrec.xgh2,...
                    satrec.xgh3,satrec.xgh4,satrec.xh2,...
                    satrec.xh3,satrec.xi2,satrec.xi3,...
                    satrec.xl2,satrec.xl3,satrec.xl4,...
                    satrec.zmol,satrec.zmos,satrec.inclo,...
                    satrec.init,ep,xincp,nodep,argpp,mp, satrec.operationmode);
                if (xincp < 0.0)
                    xincp  = -xincp;
                    nodep = nodep + pi;
                    argpp  = argpp - pi;
                end
                if ((ep < 0.0 ) || ( ep > 1.0))
                    %           fprintf(1,'# error ep %f\n', ep);
                    satrec.error = 3;
                end
            end % // if method = d

            % /* -------------------- long period periodics ------------------ */
            if (satrec.method == 'd')
                sinip =  sin(xincp);
                cosip =  cos(xincp);
                satrec.aycof = -0.5 * satrec.j3oj2 * sinip;
                % // sgp4fix for divide by zero with xinco = 180 deg
                if (abs(cosip+1.0) > 1.5e-12)
                    satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) /...
                        (1.0+cosip);
                else
                    satrec.xlcof = -0.25 * satrec.j3oj2 * sinip * (3.0 + 5.0 * cosip) /...
                        temp4;
                end
            end
            axnl = ep * cos(argpp);
            temp = 1.0 / (am * (1.0 - ep * ep));
            aynl = ep* sin(argpp) + temp * satrec.aycof;
            xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl;

            % /* --------------------- solve kepler's equation --------------- */
            u    = rem(xl - nodep, twopi);
            eo1  = u;
            tem5 = 9999.9;
            ktr = 1;
            % //   sgp4fix for kepler iteration
            % //   the following iteration needs better limits on corrections
            while (( abs(tem5) >= 1.0e-12) && (ktr <= 10) )
                sineo1 = sin(eo1);
                coseo1 = cos(eo1);
                tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
                tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
                if(abs(tem5) >= 0.95)
                    if tem5 > 0.0
                        tem5 = 0.95;
                    else
                        tem5 = -0.95;
                    end
                end
                eo1    = eo1 + tem5;
                ktr = ktr + 1;
            end

            % /* ------------- short period preliminary quantities ----------- */
            ecose = axnl*coseo1 + aynl*sineo1;
            esine = axnl*sineo1 - aynl*coseo1;
            el2   = axnl*axnl + aynl*aynl;
            pl    = am*(1.0-el2);
            if (pl < 0.0)
                %       fprintf(1,'# error pl %f\n', pl);
                satrec.error = 4;
                r = [0;0;0];
                v = [0;0;0];
            else
                rl     = am * (1.0 - ecose);
                rdotl  = sqrt(am) * esine/rl;
                rvdotl = sqrt(pl) / rl;
                betal  = sqrt(1.0 - el2);
                temp   = esine / (1.0 + betal);
                sinu   = am / rl * (sineo1 - aynl - axnl * temp);
                cosu   = am / rl * (coseo1 - axnl + aynl * temp);
                su     = atan2(sinu, cosu);
                sin2u  = (cosu + cosu) * sinu;
                cos2u  = 1.0 - 2.0 * sinu * sinu;
                temp   = 1.0 / pl;
                temp1  = 0.5 * satrec.j2 * temp;
                temp2  = temp1 * temp;

                % /* -------------- update for short period periodics ------------ */
                if (satrec.method == 'd')
                    cosisq                 = cosip * cosip;
                    satrec.con41  = 3.0*cosisq - 1.0;
                    satrec.x1mth2 = 1.0 - cosisq;
                    satrec.x7thm1 = 7.0*cosisq - 1.0;
                end
                mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) +...
                    0.5 * temp1 * satrec.x1mth2 * cos2u;
                su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
                xnode = nodep + 1.5 * temp2 * cosip * sin2u;
                xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
                mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke;
                rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +...
                    1.5 * satrec.con41) / satrec.xke;

                % /* --------------------- orientation vectors ------------------- */
                sinsu =  sin(su);
                cossu =  cos(su);
                snod  =  sin(xnode);
                cnod  =  cos(xnode);
                sini  =  sin(xinc);
                cosi  =  cos(xinc);
                xmx   = -snod * cosi;
                xmy   =  cnod * cosi;
                ux    =  xmx * sinsu + cnod * cossu;
                uy    =  xmy * sinsu + snod * cossu;
                uz    =  sini * sinsu;
                vx    =  xmx * cossu - cnod * sinsu;
                vy    =  xmy * cossu - snod * sinsu;
                vz    =  sini * cossu;

                % /* --------- position and velocity (in km and km/sec) ---------- */
                r(1) = (mrt * ux)* satrec.radiusearthkm;
                r(2) = (mrt * uy)* satrec.radiusearthkm;
                r(3) = (mrt * uz)* satrec.radiusearthkm;
                v(1) = (mvt * ux + rvdot * vx) * vkmpersec;
                v(2) = (mvt * uy + rvdot * vy) * vkmpersec;
                v(3) = (mvt * uz + rvdot * vz) * vkmpersec;
            end % // if pl > 0

            % // sgp4fix for decaying satellites
            if (mrt < 1.0)
                satrec.error = 6;
            end

        end






        function satrec = parseTle(obj, longstr1, longstr2)

            deg2rad  =   pi / 180.0;         %  0.01745329251994330;  % [deg/rad]
            xpdotp   =  1440.0 / (2.0*pi);   % 229.1831180523293;  % [rev/day]/[rad/min]

            year   = 0;
            satrec.error = 0;

            %     // set the implied decimal points since doing a formated read
            %     // fixes for bad input data values (missing, ...)
            for j = 11:16
                if (longstr1(j) == ' ')
                    longstr1(j) = '_';
                end
            end

            if (longstr1(45) ~= ' ')
                longstr1(44) = longstr1(45);
            end
            longstr1(45) = '.';

            if (longstr1(8) == ' ')
                longstr1(8) = 'U';
            end

            if (longstr1(10) == ' ')
                longstr1(10) = '.';
            end

            for j = 46:50
                if (longstr1(j) == ' ')
                    longstr1(j) = '0';
                end
            end
            if (longstr1(52) == ' ')
                longstr1(52) = '0';
            end
            if (longstr1(54) ~= ' ')
                longstr1(53) = longstr1(54);
            end
            longstr1(54) = '.';

            longstr2(26) = '.';

            for j = 27:33
                if (longstr2(j) == ' ')
                    longstr2(j) = '0';
                end
            end

            if (longstr1(63) == ' ')
                longstr1(63) = '0';
            end

            if ((length(longstr1) < 68) || (longstr1(68) == ' '))
                longstr1(68) = '0';
            end

            % parse first line
            carnumb = str2double(longstr1(1));
            satrec.satnum = str2double(longstr1(3:7));
            satrec.classification = longstr1(8);
            satrec.intldesg = longstr1(10:17);
            satrec.epochyr = str2double(longstr1(19:20));
            satrec.epochdays = str2double(longstr1(21:32));
            satrec.ndot = str2double(longstr1(34:43));
            satrec.nddot = str2double(longstr1(44:50));
            nexp = str2double(longstr1(51:52));
            satrec.bstar = str2double(longstr1(53:59));
            ibexp = str2double(longstr1(60:61));
            numb = str2double(longstr1(63));
            satrec.elnum = str2double(longstr1(65:68));

            cardnumb = str2double(longstr2(1));
            satrec.satnum = str2double(longstr2(3:7));
            satrec.inclo = str2double(longstr2(8:16));
            satrec.nodeo = str2double(longstr2(17:25));
            satrec.ecco = str2double(longstr2(26:33));
            satrec.argpo = str2double(longstr2(34:42));
            satrec.mo = str2double(longstr2(43:51));
            satrec.no_kozai = str2double(longstr2(52:63));
            satrec.revnum = str2double(longstr2(64:68));

            %     // ---- find no, ndot, nddot ----
            satrec.no_kozai   = satrec.no_kozai / xpdotp; %//* rad/min
            satrec.nddot= satrec.nddot * 10.0^nexp;
            % note the implied decimal is set when adjusting longstr1 above
            satrec.bstar= satrec.bstar * 10.0^ibexp;

            %     // ---- convert to sgp4 units ----
            %    satrec.a    = (satrec.no*tumin)^(-2/3);                % [er]
            satrec.ndot = satrec.ndot  / (xpdotp*1440.0);          % [rad/min^2]
            satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);     % [rad/min^3]

            %     // ---- find standard orbital elements ----
            satrec.inclo = satrec.inclo  * deg2rad;
            satrec.nodeo = satrec.nodeo * deg2rad;
            satrec.argpo = satrec.argpo  * deg2rad;
            satrec.mo    = satrec.mo     *deg2rad;

            %       // sgp4fix not needed here
            %    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
            %    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

            %     // ----------------------------------------------------------------
            %     // find sgp4epoch time of element set
            %     // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
            %     // and minutes from the epoch (time)
            %     // --------------------------------------------------------------

            %     // ------------- temp fix for years from 1957-2056 ----------------
            %     // ------ correct fix will occur when year is 4-digit in 2le ------
            if (satrec.epochyr < 57)
                year= satrec.epochyr + 2000;
            else
                year= satrec.epochyr + 1900;
            end
            satrec.jDayEpoch = obj.days2jDay( year, satrec.epochdays );
            sgp4epoch = satrec.jDayEpoch - 2433281.5; % days since 0 Jan 1950
            satrec = obj.initializeSgp4( satrec, sgp4epoch);
        end


        function julianDay = days2jDay (~, yr, days)
            % --------------- set up array of days in month  --------------
            for i= 1 : 12
                lmonth(i) = 31;
                if i == 2
                    lmonth(i)= 28;
                end
                if i == 4 || i == 6 || i == 9 || i == 11
                    lmonth(i)= 30;
                end
            end
            dayofyr= floor(days );
            % ----------------- find month and day of month ---------------
            if rem(yr-1900,4) == 0
                lmonth(2)= 29;
            end
            i= 1;
            inttemp= 0;
            while ( dayofyr > inttemp + lmonth(i) ) && ( i < 12 )
                inttemp= inttemp + lmonth(i);
                i= i+1;
            end
            mon= i;
            day= dayofyr - inttemp;
            % ----------------- find hours minutes and seconds ------------
            temp = (days - dayofyr )*24.0;
            hr  = fix( temp );
            temp = (temp-hr) * 60.0;
            min = fix( temp );
            sec = (temp-min) * 60.0;
            jd = 367.0 * yr  ...
                - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
                + floor( 275 * mon / 9.0 ) ...
                + day + 1721013.5;   % use - 678987.0 to go to mjd directly
            jdfrac = (sec + min * 60.0 + hr *3600.0) / 86400.0;
            % check jdfrac
            if jdfrac > 1.0
                jd = jd + floor(jdfrac);
                jdfrac = jdfrac - floor(jdfrac);
            end
            julianDay = jd + jdfrac;
        end

        function satrec = initializeSgp4(obj, satrec, epoch)

            % /* ------------------------ initialization --------------------- */
            % /* ----------- set all near earth variables to zero ------------ */
            satrec.isimp   = 0;   satrec.method = 'n'; satrec.aycof    = 0.0;
            satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
            satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
            satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
            satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
            satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
            satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
            satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot = 0.0;
            satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf  = 0.0;

            % /* ----------- set all deep space variables to zero ------------ */
            satrec.irez  = 0;   satrec.d2201 = 0.0; satrec.d2211 = 0.0;
            satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
            satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
            satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt  = 0.0;
            satrec.del1  = 0.0; satrec.del2  = 0.0; satrec.del3  = 0.0;
            satrec.didt  = 0.0; satrec.dmdt  = 0.0; satrec.dnodt = 0.0;
            satrec.domdt = 0.0; satrec.e3    = 0.0; satrec.ee2   = 0.0;
            satrec.peo   = 0.0; satrec.pgho  = 0.0; satrec.pho   = 0.0;
            satrec.pinco = 0.0; satrec.plo   = 0.0; satrec.se2   = 0.0;
            satrec.se3   = 0.0; satrec.sgh2  = 0.0; satrec.sgh3  = 0.0;
            satrec.sgh4  = 0.0; satrec.sh2   = 0.0; satrec.sh3   = 0.0;
            satrec.si2   = 0.0; satrec.si3   = 0.0; satrec.sl2   = 0.0;
            satrec.sl3   = 0.0; satrec.sl4   = 0.0; satrec.gsto  = 0.0;
            satrec.xfact = 0.0; satrec.xgh2  = 0.0; satrec.xgh3  = 0.0;
            satrec.xgh4  = 0.0; satrec.xh2   = 0.0; satrec.xh3   = 0.0;
            satrec.xi2   = 0.0; satrec.xi3   = 0.0; satrec.xl2   = 0.0;
            satrec.xl3   = 0.0; satrec.xl4   = 0.0; satrec.xlamo = 0.0;
            satrec.zmol  = 0.0; satrec.zmos  = 0.0; satrec.atime = 0.0;
            satrec.xli   = 0.0; satrec.xni   = 0.0;

            % sgp4fix add opsmode
            satrec.operationmode = "s";
            satrec.error = 0;

            % single averaged mean elements
            satrec.am = 0.0;
            satrec.em = 0.0;
            satrec.im = 0.0;
            satrec.Om = 0.0;
            satrec.mm = 0.0;
            satrec.nm = 0.0;

            %     /* -------------------- wgs-84 earth constants ----------------- */
            %     // sgp4fix identify constants and allow alternate values
            % sgp4fix switch to satrec so only one call is needed to initialize
            % constants
            %global tumin mu radiusearthkm xke j2 j3 j4 j3oj2
            satrec.mu     = 398600.5;            %// in km3 / s2
            satrec.radiusearthkm = 6378.137;     %// km
            satrec.xke    = 60.0 / sqrt(satrec.radiusearthkm*satrec.radiusearthkm*...
                satrec.radiusearthkm/satrec.mu);
            satrec.tumin  = 1.0 / satrec.xke;
            satrec.j2     =   0.00108262998905;
            satrec.j3     =  -0.00000253215306;
            satrec.j4     =  -0.00000161098761;
            satrec.j3oj2  =  satrec.j3 / satrec.j2;

            ss     = 78.0 / satrec.radiusearthkm + 1.0;
            qzms2t = ((120.0 - 78.0) / satrec.radiusearthkm)^4;

            % sgp4fix divisor for divide by zero check on inclination
            % the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
            % 1.5 e-12, so the threshold was changed to 1.5e-12 for consistancy
            temp4    =   1.5e-12;

            satrec.init = 'y';
            satrec.t    = 0.0;

            % /* -------------------- wgs-72 earth constants ----------------- */
            %     // sgp4fix identify constants and allow alternate values
            % global tumin mu radiusearthkm xke j2 j3 j4 j3oj2
            x2o3   = 2.0 / 3.0;
            %   global opsmode

            % /* ------------- calculate auxillary epoch quantities ---------- */
            eccsq  = satrec.ecco * satrec.ecco;
            omeosq = 1.0 - eccsq;
            rteosq = sqrt(omeosq);
            cosio  = cos(satrec.inclo);
            cosio2 = cosio * cosio;

            % /* ------------------ un-kozai the mean motion ----------------- */
            ak    = (satrec.xke / satrec.no_kozai)^x2o3;
            d1    = 0.75 * satrec.j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
            del   = d1 / (ak * ak);
            adel  = ak * (1.0 - del * del - del *...
                (1.0 / 3.0 + 134.0 * del * del / 81.0));
            del   = d1/(adel * adel);
            satrec.no  = satrec.no_kozai / (1.0 + del);

            ao    = (satrec.xke / satrec.no)^x2o3;
            sinio = sin(satrec.inclo);
            po    = ao * omeosq;
            con42 = 1.0 - 5.0 * cosio2;
            satrec.con41 = -con42-cosio2-cosio2;
            ainv  = 1.0 / ao;
            einv  = 1.0 / satrec.ecco;
            posq  = po * po;
            rp    = ao * (1.0 - satrec.ecco);
            satrec.method = 'n';

            % sgp4fix modern approach to finding sidereal time
            jdut1 = epoch + 2433281.5;
            twopi      = 2.0*pi;
            deg2rad    = pi/180.0;
            % ------------------------  implementation   ------------------
            tut1= ( jdut1 - 2451545.0 ) / 36525.0;
            temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1  ...
                + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;
            % 360/86400 = 1/240, to deg, to rad
            temp = rem( temp*deg2rad/240.0,twopi );
            % ------------------------ check quadrants --------------------
            if ( temp < 0.0 )
                temp = temp + twopi;
            end
            satrec.gsto = temp;

            satrec.a    = (satrec.no * satrec.tumin)^(-2.0/3.0);
            satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
            satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;

            if ((omeosq >= 0.0 ) || ( satrec.no >= 0.0))
                satrec.isimp = 0;
                if (rp < (220.0 / satrec.radiusearthkm + 1.0))
                    satrec.isimp = 1;
                end
                sfour  = ss;
                qzms24 = qzms2t;
                perige = (rp - 1.0) * satrec.radiusearthkm;

                % /* - for perigees below 156 km, s and qoms2t are altered - */
                if (perige < 156.0)
                    sfour = perige - 78.0;
                    if (perige < 98.0)
                        sfour = 20.0;
                    end
                    qzms24 = ((120.0 - sfour) / satrec.radiusearthkm)^4.0;
                    sfour  = sfour / satrec.radiusearthkm + 1.0;
                end
                pinvsq = 1.0 / posq;

                tsi  = 1.0 / (ao - sfour);
                satrec.eta  = ao * satrec.ecco * tsi;
                etasq = satrec.eta * satrec.eta;
                eeta  = satrec.ecco * satrec.eta;
                psisq = abs(1.0 - etasq);
                coef  = qzms24 * tsi^4.0;
                coef1 = coef / psisq^3.5;
                cc2   = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *...
                    (4.0 + etasq)) + 0.375 * satrec.j2 * tsi / psisq * satrec.con41 *...
                    (8.0 + 3.0 * etasq * (8.0 + etasq)));
                satrec.cc1   = satrec.bstar * cc2;
                cc3   = 0.0;
                if (satrec.ecco > 1.0e-4)
                    cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.no * sinio / satrec.ecco;
                end
                satrec.x1mth2 = 1.0 - cosio2;
                satrec.cc4    = 2.0* satrec.no * coef1 * ao * omeosq *...
                    (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *...
                    (0.5 + 2.0 * etasq) - satrec.j2 * tsi / (ao * psisq) *...
                    (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *...
                    (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *...
                    (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
                satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *...
                    (etasq + eeta) + eeta * etasq);
                cosio4 = cosio2 * cosio2;
                temp1  = 1.5 * satrec.j2 * pinvsq * satrec.no;
                temp2  = 0.5 * temp1 * satrec.j2 * pinvsq;
                temp3  = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.no;
                satrec.mdot     = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 +...
                    0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
                satrec.argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *...
                    (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +...
                    temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
                xhdot1            = -temp1 * cosio;
                satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +...
                    2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
                xpidot            =  satrec.argpdot+ satrec.nodedot;
                satrec.omgcof   = satrec.bstar * cc3 * cos(satrec.argpo);
                satrec.xmcof    = 0.0;
                if (satrec.ecco > 1.0e-4)
                    satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
                end
                satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
                satrec.t2cof   = 1.5 * satrec.cc1;

                % // sgp4fix for divide by zero with xinco = 180 deg
                if (abs(cosio+1.0) > 1.5e-12)
                    satrec.xlcof   = -0.25 * satrec.j3oj2 * sinio *...
                        (3.0 + 5.0 * cosio) / (1.0 + cosio);
                else
                    satrec.xlcof   = -0.25 * satrec.j3oj2 * sinio *...
                        (3.0 + 5.0 * cosio) / temp4;
                end
                satrec.aycof   = -0.5 * satrec.j3oj2 * sinio;
                satrec.delmo   = (1.0 + satrec.eta * cos(satrec.mo))^3;
                satrec.sinmao  = sin(satrec.mo);
                satrec.x7thm1  = 7.0 * cosio2 - 1.0;

                % /* --------------- deep space initialization ------------- */
                if ((2*pi / satrec.no) >= 225.0)
                    satrec.method = 'd';
                    satrec.isimp  = 1;
                    tc    =  0.0;
                    inclm = satrec.inclo;

                    [sinim,cosim,sinomm,cosomm,snodm,cnodm,day,satrec.e3,satrec.ee2,...
                        em,emsq,gam,satrec.peo,satrec.pgho,satrec.pho,satrec.pinco,...
                        satrec.plo,rtemsq,satrec.se2,satrec.se3,satrec.sgh2,...
                        satrec.sgh3,satrec.sgh4,satrec.sh2,satrec.sh3,satrec.si2,...
                        satrec.si3,satrec.sl2,satrec.sl3,satrec.sl4,s1,s2,s3,s4,s5,...
                        s6,s7,ss1,ss2,ss3,ss4,ss5,ss6,ss7,sz1,sz2,sz3,sz11,sz12,...
                        sz13,sz21,sz22,sz23,sz31,sz32,sz33,satrec.xgh2,satrec.xgh3,...
                        satrec.xgh4,satrec.xh2,satrec.xh3,satrec.xi2,satrec.xi3,...
                        satrec.xl2,satrec.xl3,satrec.xl4,nm,z1,z2,z3,z11,z12,z13,...
                        z21,z22,z23,z31,z32,z33,satrec.zmol,satrec.zmos] ...
                        = obj.dscom(epoch,satrec.ecco,satrec.argpo,tc,satrec.inclo,...
                        satrec.nodeo,satrec.no);

                    [satrec.ecco,satrec.inclo,satrec.nodeo,satrec.argpo,satrec.mo]...
                        = obj.dpper(satrec.e3,satrec.ee2,satrec.peo,satrec.pgho,...
                        satrec.pho,satrec.pinco,satrec.plo,satrec.se2,satrec.se3,...
                        satrec.sgh2,satrec.sgh3,satrec.sgh4,satrec.sh2,satrec.sh3,...
                        satrec.si2,satrec.si3,satrec.sl2,satrec.sl3,satrec.sl4,...
                        satrec.t,satrec.xgh2,satrec.xgh3,satrec.xgh4,satrec.xh2,...
                        satrec.xh3,satrec.xi2,satrec.xi3,satrec.xl2,satrec.xl3,...
                        satrec.xl4,satrec.zmol,satrec.zmos,inclm,satrec.init,...
                        satrec.ecco,satrec.inclo,satrec.nodeo,satrec.argpo,satrec.mo, ...
                        satrec.operationmode);

                    argpm  = 0.0;
                    nodem  = 0.0;
                    mm     = 0.0;

                    [em,argpm,inclm,mm,nm,nodem,satrec.irez,satrec.atime,...
                        satrec.d2201,satrec.d2211,satrec.d3210,satrec.d3222,...
                        satrec.d4410,satrec.d4422,satrec.d5220,satrec.d5232,...
                        satrec.d5421,satrec.d5433,satrec.dedt,satrec.didt,...
                        satrec.dmdt,dndt,satrec.dnodt,satrec.domdt,satrec.del1,...
                        satrec.del2,satrec.del3,...
                        satrec.xfact,satrec.xlamo,satrec.xli,satrec.xni] ...
                        = obj.dsinit(...
                        satrec.xke, cosim,emsq,satrec.argpo,s1,s2,s3,s4,s5,sinim,ss1,ss2,ss3,...
                        ss4,ss5,sz1,sz3,sz11,sz13,sz21,sz23,sz31,sz33,satrec.t,tc,...
                        satrec.gsto,satrec.mo,satrec.mdot,satrec.no,satrec.nodeo,...
                        satrec.nodedot,xpidot,z1,z3,z11,z13,z21,z23,z31,z33,em,...
                        argpm,inclm,mm,nm,nodem,satrec.ecco,eccsq);
                end

                % /* ----------- set variables if not deep space ----------- */
                if (satrec.isimp ~= 1)
                    cc1sq          = satrec.cc1 * satrec.cc1;
                    satrec.d2    = 4.0 * ao * tsi * cc1sq;
                    temp           = satrec.d2 * tsi * satrec.cc1 / 3.0;
                    satrec.d3    = (17.0 * ao + sfour) * temp;
                    satrec.d4    = 0.5 * temp * ao * tsi *...
                        (221.0 * ao + 31.0 * sfour) * satrec.cc1;
                    satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
                    satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *...
                        (12.0 * satrec.d2 + 10.0 * cc1sq));
                    satrec.t5cof = 0.2 * (3.0 * satrec.d4 +...
                        12.0 * satrec.cc1 * satrec.d3 +...
                        6.0 * satrec.d2 * satrec.d2 +...
                        15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
                end
            end % // if omeosq = 0 ...

            % /* finally propogate to zero epoch to initialise all others. */
            % sgp4fix take out check to let satellites process until they are actually below earth surface
            %if(satrec.error == 0)
            % [satrec, r, v] = sgp4(satrec, 0.0);
            %end

            satrec.init = 'n';

        end


        function [  ep,     inclp,  nodep, argpp,  mp]...
                = dpper(...
                ~, e3,     ee2,    peo,    pgho,   pho,    pinco,  plo,    se2,...
                se3,    sgh2,   sgh3,   sgh4,   sh2,    sh3,    si2,    si3,...
                sl2,    sl3,    sl4,    t,      xgh2,   xgh3,   xgh4,   xh2,...
                xh3,    xi2,    xi3,    xl2,    xl3,    xl4,    zmol,...
                zmos,   inclo,  init,   ep,     inclp,  nodep, argpp,  mp, opsmode)

            % change to variable passed in
            % global opsmode

            % /* --------------------- local variables ------------------------ */
            twopi = 2.0 * pi;

            % /* ---------------------- constants ----------------------------- */
            zns   = 1.19459e-5;
            zes   = 0.01675;
            znl   = 1.5835218e-4;
            zel   = 0.05490;

            % /* --------------- calculate time varying periodics ----------- */
            zm    = zmos + zns * t;
            % // be sure that the initial call has time set to zero
            if (init == 'y')
                zm = zmos;
            end
            zf    = zm + 2.0 * zes * sin(zm);
            sinzf = sin(zf);
            f2    =  0.5 * sinzf * sinzf - 0.25;
            f3    = -0.5 * sinzf * cos(zf);
            ses   = se2* f2 + se3 * f3;
            sis   = si2 * f2 + si3 * f3;
            sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
            sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
            shs   = sh2 * f2 + sh3 * f3;
            zm    = zmol + znl * t;
            if (init == 'y')
                zm = zmol;
            end
            zf    = zm + 2.0 * zel * sin(zm);
            sinzf = sin(zf);
            f2    =  0.5 * sinzf * sinzf - 0.25;
            f3    = -0.5 * sinzf * cos(zf);
            sel   = ee2 * f2 + e3 * f3;
            sil   = xi2 * f2 + xi3 * f3;
            sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
            sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
            shll  = xh2 * f2 + xh3 * f3;
            pe    = ses + sel;
            pinc  = sis + sil;
            pl    = sls + sll;
            pgh   = sghs + sghl;
            ph    = shs + shll;

            if (init == 'n')
                %   //  0.2 rad = 11.45916 deg
                pe    = pe - peo;
                pinc  = pinc - pinco;
                pl    = pl - plo;
                pgh   = pgh - pgho;
                ph    = ph - pho;
                inclp = inclp + pinc;
                ep    = ep + pe;
                sinip = sin(inclp);
                cosip = cos(inclp);

                % /* ----------------- apply periodics directly ------------ */
                %   //  sgp4fix for lyddane choice
                %   //  strn3 used original inclination - this is technically feasible
                %   //  gsfc used perturbed inclination - also technically feasible
                %   //  probably best to readjust the 0.2 limit value and limit discontinuity
                %   //  use next line for original strn3 approach and original inclination
                %   //  if (inclo >= 0.2)
                %   //  use next line for gsfc version and perturbed inclination
                if (inclp >= 0.2)
                    ph     = ph / sinip;
                    pgh    = pgh - cosip * ph;
                    argpp  = argpp + pgh;
                    nodep  = nodep + ph;
                    mp     = mp + pl;
                else
                    % /* ---- apply periodics with lyddane modification ---- */
                    sinop  = sin(nodep);
                    cosop  = cos(nodep);
                    alfdp  = sinip * sinop;
                    betdp  = sinip * cosop;
                    dalf   =  ph * cosop + pinc * cosip * sinop;
                    dbet   = -ph * sinop + pinc * cosip * cosop;
                    alfdp  = alfdp + dalf;
                    betdp  = betdp + dbet;
                    nodep  = rem(nodep, twopi);
                    % sgp4fix for afspc written intrinsic functions
                    % nodep used without a trigonometric function ahead
                    if ((nodep < 0.0) && (opsmode == 'a'))
                        nodep = nodep + twopi;
                    end
                    xls    = mp + argpp + cosip * nodep;
                    dls    = pl + pgh - pinc * nodep * sinip;
                    xls    = xls + dls;
                    xnoh   = nodep;
                    nodep  = atan2(alfdp, betdp);
                    % sgp4fix for afspc written intrinsic functions
                    % nodep used without a trigonometric function ahead
                    if ((nodep < 0.0) && (opsmode == 'a'))
                        nodep = nodep + twopi;
                    end
                    if (abs(xnoh - nodep) > pi)
                        if (nodep < xnoh)
                            nodep = nodep + twopi;
                        else
                            nodep = nodep - twopi;
                        end
                    end
                    mp    = mp + pl;
                    argpp = xls - mp - cosip * nodep;
                end
            end % // if init == 'n'

        end

        function [ sinim,cosim,sinomm,cosomm,snodm,cnodm,day,e3,ee2,em,emsq,gam,...
                peo,pgho,pho,pinco,plo,rtemsq,se2,se3,sgh2,sgh3,sgh4,sh2,sh3,si2,...
                si3,sl2,sl3,sl4,s1,s2,s3,s4,s5,s6,s7,ss1,ss2,ss3,ss4,ss5,ss6,ss7,...
                sz1,sz2,sz3,sz11,sz12,sz13,sz21,sz22,sz23,sz31,sz32,sz33,xgh2,xgh3,...
                xgh4,xh2,xh3,xi2,xi3,xl2,xl3,xl4,nm,z1,z2,z3,z11,z12,z13,z21,z22,...
                z23,z31,z32,z33,zmol,zmos]...
                = dscom (~, epoch, ep, argpp, tc, inclp, nodep, np)

            % /* -------------------------- constants ------------------------- */
            zes     =  0.01675;
            zel     =  0.05490;
            c1ss    =  2.9864797e-6;
            c1l     =  4.7968065e-7;
            zsinis  =  0.39785416;
            zcosis  =  0.91744867;
            zcosgs  =  0.1945905;
            zsings  = -0.98088458;
            twopi   =  2.0 * pi;

            % /* --------------------- local variables ------------------------ */
            nm     = np;
            em     = ep;
            snodm  = sin(nodep);
            cnodm  = cos(nodep);
            sinomm = sin(argpp);
            cosomm = cos(argpp);
            sinim  = sin(inclp);
            cosim  = cos(inclp);
            emsq   = em * em;
            betasq = 1.0 - emsq;
            rtemsq = sqrt(betasq);

            % /* ----------------- initialize lunar solar terms --------------- */
            peo    = 0.0;
            pinco  = 0.0;
            plo    = 0.0;
            pgho   = 0.0;
            pho    = 0.0;
            day    = epoch + 18261.5 + tc / 1440.0;
            xnodce = rem(4.5236020 - 9.2422029e-4 * day, twopi);
            stem   = sin(xnodce);
            ctem   = cos(xnodce);
            zcosil = 0.91375164 - 0.03568096 * ctem;
            zsinil = sqrt(1.0 - zcosil * zcosil);
            zsinhl = 0.089683511 * stem / zsinil;
            zcoshl = sqrt(1.0 - zsinhl * zsinhl);
            gam    = 5.8351514 + 0.0019443680 * day;
            zx     = 0.39785416 * stem / zsinil;
            zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
            zx     = atan2(zx, zy);
            zx     = gam + zx - xnodce;
            zcosgl = cos(zx);
            zsingl = sin(zx);

            % /* ------------------------- do solar terms --------------------- */
            zcosg = zcosgs;
            zsing = zsings;
            zcosi = zcosis;
            zsini = zsinis;
            zcosh = cnodm;
            zsinh = snodm;
            cc    = c1ss;
            xnoi  = 1.0 / nm;

            for (lsflg = 1:2)
                a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
                a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
                a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
                a8  =   zsing * zsini;
                a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
                a10 =   zcosg * zsini;
                a2  =   cosim * a7 + sinim * a8;
                a4  =   cosim * a9 + sinim * a10;
                a5  =  -sinim * a7 + cosim * a8;
                a6  =  -sinim * a9 + cosim * a10;

                x1  =  a1 * cosomm + a2 * sinomm;
                x2  =  a3 * cosomm + a4 * sinomm;
                x3  = -a1 * sinomm + a2 * cosomm;
                x4  = -a3 * sinomm + a4 * cosomm;
                x5  =  a5 * sinomm;
                x6  =  a6 * sinomm;
                x7  =  a5 * cosomm;
                x8  =  a6 * cosomm;

                z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
                z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
                z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
                z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq;
                z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq;
                z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq;
                z11 = -6.0 * a1 * a5 + emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
                z12 = -6.0 *  (a1 * a6 + a3 * a5) + emsq *...
                    (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
                z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
                z21 =  6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
                z22 =  6.0 *  (a4 * a5 + a2 * a6) + emsq *...
                    (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
                z23 =  6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
                z1  = z1 + z1 + betasq * z31;
                z2  = z2 + z2 + betasq * z32;
                z3  = z3 + z3 + betasq * z33;
                s3  = cc * xnoi;
                s2  = -0.5 * s3 / rtemsq;
                s4  = s3 * rtemsq;
                s1  = -15.0 * em * s4;
                s5  = x1 * x3 + x2 * x4;
                s6  = x2 * x3 + x1 * x4;
                s7  = x2 * x4 - x1 * x3;

                % /* ----------------------- do lunar terms ------------------- */
                if (lsflg == 1)
                    ss1   = s1;
                    ss2   = s2;
                    ss3   = s3;
                    ss4   = s4;
                    ss5   = s5;
                    ss6   = s6;
                    ss7   = s7;
                    sz1   = z1;
                    sz2   = z2;
                    sz3   = z3;
                    sz11  = z11;
                    sz12  = z12;
                    sz13  = z13;
                    sz21  = z21;
                    sz22  = z22;
                    sz23  = z23;
                    sz31  = z31;
                    sz32  = z32;
                    sz33  = z33;
                    zcosg = zcosgl;
                    zsing = zsingl;
                    zcosi = zcosil;
                    zsini = zsinil;
                    zcosh = zcoshl * cnodm + zsinhl * snodm;
                    zsinh = snodm * zcoshl - cnodm * zsinhl;
                    cc    = c1l;
                end
            end

            zmol = rem(4.7199672 + 0.22997150  * day - gam, twopi);
            zmos = rem(6.2565837 + 0.017201977 * day, twopi);

            % /* ------------------------ do solar terms ---------------------- */
            se2  =   2.0 * ss1 * ss6;
            se3  =   2.0 * ss1 * ss7;
            si2  =   2.0 * ss2 * sz12;
            si3  =   2.0 * ss2 * (sz13 - sz11);
            sl2  =  -2.0 * ss3 * sz2;
            sl3  =  -2.0 * ss3 * (sz3 - sz1);
            sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
            sgh2 =   2.0 * ss4 * sz32;
            sgh3 =   2.0 * ss4 * (sz33 - sz31);
            sgh4 = -18.0 * ss4 * zes;
            sh2  =  -2.0 * ss2 * sz22;
            sh3  =  -2.0 * ss2 * (sz23 - sz21);

            % /* ------------------------ do lunar terms ---------------------- */
            ee2  =   2.0 * s1 * s6;
            e3   =   2.0 * s1 * s7;
            xi2  =   2.0 * s2 * z12;
            xi3  =   2.0 * s2 * (z13 - z11);
            xl2  =  -2.0 * s3 * z2;
            xl3  =  -2.0 * s3 * (z3 - z1);
            xl4  =  -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
            xgh2 =   2.0 * s4 * z32;
            xgh3 =   2.0 * s4 * (z33 - z31);
            xgh4 = -18.0 * s4 * zel;
            xh2  =  -2.0 * s2 * z22;
            xh3  =  -2.0 * s2 * (z23 - z21);

        end


        function [  em,     argpm,  inclm,  mm,     nm,     nodem, irez,...
                atime,  d2201,  d2211,  d3210,  d3222,  d4410,  d4422,...
                d5220,  d5232,  d5421,  d5433,  dedt,   didt,   dmdt,...
                dndt,   dnodt,  domdt,  del1,   del2,   del3,   xfact,...
                xlamo,  xli,    xni]...
                = dsinit( ...
                ~, xke,    cosim,  emsq,   argpo,  s1,     s2,     s3,     s4,...
                s5,     sinim,  ss1,    ss2,    ss3,    ss4,    ss5,...
                sz1,    sz3,    sz11,   sz13,   sz21,   sz23,   sz31,...
                sz33,   t,      tc,     gsto,   mo,     mdot,   no,...
                nodeo,  nodedot,xpidot, z1,     z3,     z11,...
                z13,    z21,    z23,    z31,    z33,    em,     argpm,...
                inclm,  mm,     nm,     nodem,  ecco,   eccsq)

            % /* --------------------- local variables ------------------------ */
            twopi = 2.0 * pi;
            aonv  = 0.0;
            q22    = 1.7891679e-6;
            q31    = 2.1460748e-6;
            q33    = 2.2123015e-7;
            root22 = 1.7891679e-6;
            root44 = 7.3636953e-9;
            root54 = 2.1765803e-9;
            rptim  = 4.37526908801129966e-3;
            root32 = 3.7393792e-7;
            root52 = 1.1428639e-7;
            x2o3   = 2.0 / 3.0;
            znl    = 1.5835218e-4;
            zns    = 1.19459e-5;

            %     // sgp4fix identify constants and allow alternate values
            % sgp4fix no longer needed, pass xke in
            %global tumin mu radiusearthkm xke j2 j3 j4 j3oj2

            % /* -------------------- deep space initialization ------------ */
            irez = 0;
            if ((nm < 0.0052359877) && (nm > 0.0034906585))
                irez = 1;
            end
            if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5))
                irez = 2;
            end
            d2201 = 0;
            d2211 = 0;
            d3210 = 0;
            d3222 = 0;
            d4410 = 0;
            d4422 = 0;
            d5220 = 0;
            d5232 = 0;
            d5421 = 0;
            d5433 = 0;
            del1  = 0;
            del2  = 0;
            del3  = 0;
            atime = 0;
            xfact = 0;
            xlamo = 0;
            xli   = 0;
            xni   = 0;

            % /* ------------------------ do solar terms ------------------- */
            ses  =  ss1 * zns * ss5;
            sis  =  ss2 * zns * (sz11 + sz13);
            sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
            sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
            shs  = -zns * ss2 * (sz21 + sz23);
            %   // sgp4fix for 180 deg incl
            if ((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2))
                shs = 0.0;
            end
            if (sinim ~= 0.0)
                shs = shs / sinim;
            end
            sgs  = sghs - cosim * shs;

            % /* ------------------------- do lunar terms ------------------ */
            dedt = ses + s1 * znl * s5;
            didt = sis + s2 * znl * (z11 + z13);
            dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
            sghl = s4 * znl * (z31 + z33 - 6.0);
            shll = -znl * s2 * (z21 + z23);
            %   // sgp4fix for 180 deg incl
            if ((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2))
                shll = 0.0;
            end
            domdt = sgs + sghl;
            dnodt = shs;
            if (sinim ~= 0.0)
                domdt = domdt - cosim / sinim * shll;
                dnodt = dnodt + shll / sinim;
            end

            % /* ----------- calculate deep space resonance effects -------- */
            dndt   = 0.0;
            theta  = rem(gsto + tc * rptim, twopi);
            em     = em + dedt * t;
            inclm  = inclm + didt * t;
            argpm  = argpm + domdt * t;
            nodem  = nodem + dnodt * t;
            mm     = mm + dmdt * t;
            % //   sgp4fix for negative inclinations
            % //   the following if statement should be commented out
            % //if (inclm < 0.0)
            % //  {
            % //    inclm  = -inclm;
            % //    argpm  = argpm - pi;
            % //    nodem = nodem + pi;
            % //  }

            %  /* - update resonances : numerical (euler-maclaurin) integration - */
            %  /* ------------------------- epoch restart ----------------------  */
            %  //   sgp4fix for propagator problems
            %  //   the following integration works for negative time steps and periods
            %  //   the specific changes are unknown because the original code was so convoluted

            % /* -------------- initialize the resonance terms ------------- */
            if (irez ~= 0)
                aonv = (nm / xke)^x2o3;

                % /* ---------- geopotential resonance for 12 hour orbits ------ */
                if (irez == 2)
                    cosisq = cosim * cosim;
                    emo    = em;
                    em     = ecco;
                    emsqo  = emsq;
                    emsq   = eccsq;
                    eoc    = em * emsq;
                    g201   = -0.306 - (em - 0.64) * 0.440;

                    if (em <= 0.65)
                        g211 =    3.616  -  13.2470 * em +  16.2900 * emsq;
                        g310 =  -19.302  + 117.3900 * em - 228.4190 * emsq +  156.5910 * eoc;
                        g322 =  -18.9068 + 109.7927 * em - 214.6334 * emsq +  146.5816 * eoc;
                        g410 =  -41.122  + 242.6940 * em - 471.0940 * emsq +  313.9530 * eoc;
                        g422 = -146.407  + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
                        g520 = -532.114  + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
                    else
                        g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
                        g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
                        g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
                        g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
                        g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
                        if (em > 0.715)
                            g520 =-5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
                        else
                            g520 = 1464.74 -  4664.75 * em +  3763.64 * emsq;
                        end
                    end
                    if (em < 0.7)
                        g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
                        g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
                        g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc;
                    else
                        g533 =-37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
                        g521 =-51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
                        g532 =-40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
                    end

                    sini2=  sinim * sinim;
                    f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
                    f221 =  1.5 * sini2;
                    f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
                    f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
                    f441 = 35.0 * sini2 * f220;
                    f442 = 39.3750 * sini2 * sini2;
                    f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) +...
                        0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
                    f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +...
                        10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
                    f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *...
                        (-12.0 + 8.0 * cosim + 10.0 * cosisq));
                    f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *...
                        (12.0 + 8.0 * cosim - 10.0 * cosisq));
                    xno2  =  nm * nm;
                    ainv2 =  aonv * aonv;
                    temp1 =  3.0 * xno2 * ainv2;
                    temp  =  temp1 * root22;
                    d2201 =  temp * f220 * g201;
                    d2211 =  temp * f221 * g211;
                    temp1 =  temp1 * aonv;
                    temp  =  temp1 * root32;
                    d3210 =  temp * f321 * g310;
                    d3222 =  temp * f322 * g322;
                    temp1 =  temp1 * aonv;
                    temp  =  2.0 * temp1 * root44;
                    d4410 =  temp * f441 * g410;
                    d4422 =  temp * f442 * g422;
                    temp1 =  temp1 * aonv;
                    temp  =  temp1 * root52;
                    d5220 =  temp * f522 * g520;
                    d5232 =  temp * f523 * g532;
                    temp  =  2.0 * temp1 * root54;
                    d5421 =  temp * f542 * g521;
                    d5433 =  temp * f543 * g533;
                    xlamo =  rem(mo + nodeo + nodeo-theta - theta, twopi);
                    xfact =  mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
                    em    = emo;
                    emsq  = emsqo;
                end

                % /* ---------------- synchronous resonance terms -------------- */
                if (irez == 1)
                    g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
                    g310  = 1.0 + 2.0 * emsq;
                    g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
                    f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim);
                    f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
                    f330  = 1.0 + cosim;
                    f330  = 1.875 * f330 * f330 * f330;
                    del1  = 3.0 * nm * nm * aonv * aonv;
                    del2  = 2.0 * del1 * f220 * g200 * q22;
                    del3  = 3.0 * del1 * f330 * g300 * q33 * aonv;
                    del1  = del1 * f311 * g310 * q31 * aonv;
                    xlamo = rem(mo + nodeo + argpo - theta, twopi);
                    xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;
                end

                % /* ------------ for sgp4, initialize the integrator ---------- */
                xli   = xlamo;
                xni   = no;
                atime = 0.0;
                nm    = no + dndt;
            end

        end

        function [  atime,  em,     argpm,  inclm,  xli,    mm,     xni,...
                nodem, dndt,   nm]...
                = dspace(~, ...
                d2201,  d2211,  d3210,  d3222,  d4410,  d4422,  d5220,...
                d5232,  d5421,  d5433,  dedt,   del1,   del2,   del3,...
                didt,   dmdt,   dnodt,  domdt,  irez,   argpo,  argpdot,...
                t,      tc,     gsto,   xfact,  xlamo,  no,     atime,...
                em,     argpm,  inclm,  xli,    mm,     xni,    nodem,...
                nm)

            twopi = 2.0 * pi;

            fasx2 = 0.13130908;
            fasx4 = 2.8843198;
            fasx6 = 0.37448087;
            g22   = 5.7686396;
            g32   = 0.95240898;
            g44   = 1.8014998;
            g52   = 1.0508330;
            g54   = 4.4108898;
            rptim = 4.37526908801129966e-3;
            stepp =    720.0;
            stepn =   -720.0;
            step2 = 259200.0;

            % /* ----------- calculate deep space resonance effects ----------- */
            dndt   = 0.0;
            theta  = rem(gsto + tc * rptim, twopi);
            em     = em + dedt * t;

            inclm  = inclm + didt * t;
            argpm  = argpm + domdt * t;
            nodem  = nodem + dnodt * t;
            mm     = mm + dmdt * t;

            % //   sgp4fix for negative inclinations
            % //   the following if statement should be commented out
            % //  if (inclm < 0.0)
            % // {
            % //    inclm  = -inclm;
            % //    argpm  = argpm - pi;
            % //    nodem = nodem + pi;
            % //  }

            % /* - update resonances : numerical (euler-maclaurin) integration - */
            % /* ------------------------- epoch restart ----------------------  */

            % //   sgp4fix for propagator problems
            % //   the following integration works for negative time steps and periods
            % //   the specific changes are unknown because the original code was so convoluted

            % // sgp4fix take out atime = 0.0 and fix for faster operation
            ft    = 0.0;

            if (irez ~= 0)
                % sgp4fix streamline check
                if ((atime == 0.0) || (t * atime <= 0.0) || ...
                        (abs(t) < abs(atime)) )
                    atime  = 0.0;
                    xni    = no;
                    xli    = xlamo;
                end
                % sgp4fix move check outside loop
                if (t >= 0.0)
                    delt = stepp;
                else
                    delt = stepn;
                end

                iretn = 381; %// added for do loop
                iret  =   0; %// added for loop
                while (iretn == 381)
                    % /* ------------------- dot terms calculated ------------- */
                    % /* ----------- near - synchronous resonance terms ------- */
                    if (irez ~= 2)
                        xndt  = del1 * sin(xli - fasx2) + del2 * sin(2.0 * (xli - fasx4)) +...
                            del3 * sin(3.0 * (xli - fasx6));
                        xldot = xni + xfact;
                        xnddt = del1 * cos(xli - fasx2) +...
                            2.0 * del2 * cos(2.0 * (xli - fasx4)) +...
                            3.0 * del3 * cos(3.0 * (xli - fasx6));
                        xnddt = xnddt * xldot;
                    else
                        % /* --------- near - half-day resonance terms -------- */
                        xomi  = argpo + argpdot * atime;
                        x2omi = xomi + xomi;
                        x2li  = xli + xli;
                        xndt  = d2201 * sin(x2omi + xli - g22) + d2211 * sin(xli - g22) +...
                            d3210 * sin(xomi + xli - g32)  + d3222 * sin(-xomi + xli - g32)+...
                            d4410 * sin(x2omi + x2li - g44)+ d4422 * sin(x2li - g44) +...
                            d5220 * sin(xomi + xli - g52)  + d5232 * sin(-xomi + xli - g52)+...
                            d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54);
                        xldot = xni + xfact;
                        xnddt = d2201 * cos(x2omi + xli - g22) + d2211 * cos(xli - g22) +...
                            d3210 * cos(xomi + xli - g32) + d3222 * cos(-xomi + xli - g32) +...
                            d5220 * cos(xomi + xli - g52) + d5232 * cos(-xomi + xli - g52) +...
                            2.0 * (d4410 * cos(x2omi + x2li - g44) +...
                            d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +...
                            d5433 * cos(-xomi + x2li - g54));
                        xnddt = xnddt * xldot;
                    end

                    % /* ----------------------- integrator ------------------- */
                    % sgp4fix move end checks to end of routine
                    if (abs(t - atime) >= stepp)
                        iret  = 0;
                        iretn = 381;
                    else
                        ft    = t - atime;
                        iretn = 0;
                    end

                    if (iretn == 381)
                        xli   = xli + xldot * delt + xndt * step2;
                        xni   = xni + xndt * delt + xnddt * step2;
                        atime = atime + delt;
                    end
                end %  // while iretn = 381

                nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
                xl = xli + xldot * ft + xndt * ft * ft * 0.5;
                if (irez ~= 1)
                    mm   = xl - 2.0 * nodem + 2.0 * theta;
                    dndt = nm - no;
                else
                    mm   = xl - nodem - argpm+ theta;
                    dndt = nm - no;
                end

                nm = no + dndt;
            end
        end

        function [prec,psia,wa,ea,xa] = precess (~, ttt )

            % " to rad
            convrt = pi / (180.0*3600.0);
            ttt2= ttt * ttt;
            ttt3= ttt2 * ttt;

            prec = eye(3);

            %     fprintf(1,'80prec %15.9f  \n',ttt);
            psia =             5038.7784*ttt - 1.07259*ttt2 - 0.001147*ttt3; % "
            wa   = 84381.448                 + 0.05127*ttt2 - 0.007726*ttt3;
            ea   = 84381.448 -   46.8150*ttt - 0.00059*ttt2 + 0.001813*ttt3;
            xa   =               10.5526*ttt - 2.38064*ttt2 - 0.001125*ttt3;

            zeta =             2306.2181*ttt + 0.30188*ttt2 + 0.017998*ttt3; % "
            theta=             2004.3109*ttt - 0.42665*ttt2 - 0.041833*ttt3;
            z    =             2306.2181*ttt + 1.09468*ttt2 + 0.018203*ttt3;
            % ------------------ iau 06 precession angles -------------------


            % convert units to rad
            psia = psia  * convrt; % rad
            wa   = wa    * convrt;
            ea   = ea    * convrt;
            xa   = xa    * convrt;

            zeta = zeta  * convrt;
            theta= theta * convrt;
            z    = z     * convrt;

            coszeta  = cos(zeta);
            sinzeta  = sin(zeta);
            costheta = cos(theta);
            sintheta = sin(theta);
            cosz     = cos(z);
            sinz     = sin(z);

            % ----------------- form matrix  mod to j2000 -----------------
            prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz;
            prec(1,2) =  coszeta * costheta * sinz + sinzeta * cosz;
            prec(1,3) =  coszeta * sintheta;
            prec(2,1) = -sinzeta * costheta * cosz - coszeta * sinz;
            prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz;
            prec(2,3) = -sinzeta * sintheta;
            prec(3,1) = -sintheta * cosz;
            prec(3,2) = -sintheta * sinz;
            prec(3,3) =  costheta;
        end

        function [deltapsi, trueeps, meaneps, omega, nut] = nutation(obj, ttt, ddpsi, ddeps)

            deg2rad = pi/180.0;
            iar80 = obj.iar80Ar;
            rar80 = obj.rar80Ar;

            % ---- determine coefficients for iau 1980 nutation theory ----
            ttt2= ttt*ttt;
            ttt3= ttt2*ttt;

            meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
            meaneps = rem( meaneps/3600.0, 360.0 );
            meaneps = meaneps * deg2rad;

            l = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) / 3600.0 + 134.96298139;
            l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333;
            f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028;
            d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306;
            omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222;

            % ---- convert units to rad
            l    = rem( l,360.0  )     * deg2rad; % rad
            l1   = rem( l1,360.0  )    * deg2rad;
            f    = rem( f,360.0  )     * deg2rad;
            d    = rem( d,360.0  )     * deg2rad;
            omega= rem( omega,360.0  ) * deg2rad;

            deltapsi= 0.0;
            deltaeps= 0.0;
            for i= 106:-1: 1
                tempval= iar80(i,1)*l + iar80(i,2)*l1 + iar80(i,3)*f + ...
                    iar80(i,4)*d + iar80(i,5)*omega;
                deltapsi= deltapsi + (rar80(i,1)+rar80(i,2)*ttt) * sin( tempval );
                deltaeps= deltaeps + (rar80(i,3)+rar80(i,4)*ttt) * cos( tempval );
            end

            % --------------- find nutation parameters --------------------
            deltapsi = rem( deltapsi + ddpsi, 2.0 * pi );
            deltaeps = rem( deltaeps + ddeps, 2.0 * pi );
            trueeps  = meaneps + deltaeps;

            cospsi  = cos(deltapsi);
            sinpsi  = sin(deltapsi);
            coseps  = cos(meaneps);
            sineps  = sin(meaneps);
            costrueeps = cos(trueeps);
            sintrueeps = sin(trueeps);

            nut(1,1) =  cospsi;
            nut(1,2) =  costrueeps * sinpsi;
            nut(1,3) =  sintrueeps * sinpsi;
            nut(2,1) = -coseps * sinpsi;
            nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps;
            nut(2,3) =  sintrueeps * coseps * cospsi - sineps * costrueeps;
            nut(3,1) = -sineps * sinpsi;
            nut(3,2) =  costrueeps * sineps * cospsi - sintrueeps * coseps;
            nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps;
        end

        function [reci, veci] = temeToEci(obj, rteme, vteme, jd)
            ttt = (jd - 2451545.0  )/ 36525.0;
            [ddpsi, ddeps] = obj.getPsiEps(jd);
            [prec,psia, wa, ea, xa] = obj.precess( ttt );
            [deltapsi, trueeps, meaneps, omega, nut] = obj.nutation(ttt, ddpsi, ddeps );
            % ------------------------ find eqeg ----------------------
            % rotate teme through just geometric terms
            eqeg = deltapsi* cos(meaneps);
            eqeg = rem (eqeg, 2.0*pi);
            eqe(1,1) =  cos(eqeg);
            eqe(1,2) =  sin(eqeg);
            eqe(1,3) =  0.0;
            eqe(2,1) = -sin(eqeg);
            eqe(2,2) =  cos(eqeg);
            eqe(2,3) =  0.0;
            eqe(3,1) =  0.0;
            eqe(3,2) =  0.0;
            eqe(3,3) =  1.0;
            tm = prec * nut * eqe';
            reci = tm * rteme';
            veci = tm * vteme';
        end
        function [ddpsi, ddeps] = getPsiEps(obj, jd)
            mjd = jd - 2400000.5; % modify julian Day
            eopData = obj.getRecentEop()';
            firstDay = eopData(1,4);
            lastDay = eopData(end,4);
            conv = pi / (180*3600);
            if mjd < firstDay
                ddpsi = eopData(1,9) * conv;
                ddeps = eopData(1,10) * conv;
            elseif mjd > lastDay
                ddpsi = eopData(end,9) * conv;
                ddeps = eopData(end,10) * conv;
            else
                ddpsi = eopData(floor(mjd-firstDay) + 1,9) * conv;
                ddeps = eopData(floor(mjd-firstDay) + 1,10) * conv;
            end

        end

        function nut80 = getNut80(~)
            nut80 = [0	0	0	0	1	-171996	-174.2	92025	8.9	1
                0	0	2	-2	2	-13187	-1.6	5736	-3.1	9
                0	0	2	0	2	-2274	-0.2	977	-0.5	31
                0	0	0	0	2	2062	0.2	-895	0.5	2
                0	1	0	0	0	1426	-3.4	54	-0.1	10
                1	0	0	0	0	712	0.1	-7	0	32
                0	1	2	-2	2	-517	1.2	224	-0.6	11
                0	0	2	0	1	-386	-0.4	200	0	33
                1	0	2	0	2	-301	0	129	-0.1	34
                0	-1	2	-2	2	217	-0.5	-95	0.3	12
                1	0	0	-2	0	-158	0	-1	0	35
                0	0	2	-2	1	129	0.1	-70	0	13
                -1	0	2	0	2	123	0	-53	0	36
                1	0	0	0	1	63	0.1	-33	0	38
                0	0	0	2	0	63	0	-2	0	37
                -1	0	2	2	2	-59	0	26	0	40
                -1	0	0	0	1	-58	-0.1	32	0	39
                1	0	2	0	1	-51	0	27	0	41
                2	0	0	-2	0	48	0	1	0	14
                -2	0	2	0	1	46	0	-24	0	3
                0	0	2	2	2	-38	0	16	0	42
                2	0	2	0	2	-31	0	13	0	45
                2	0	0	0	0	29	0	-1	0	43
                1	0	2	-2	2	29	0	-12	0	44
                0	0	2	0	0	26	0	-1	0	46
                0	0	2	-2	0	-22	0	0	0	15
                -1	0	2	0	1	21	0	-10	0	47
                0	2	0	0	0	17	-0.100000000000000	0	0	16
                0	2	2	-2	2	-16	0.100000000000000	7	0	18
                -1	0	0	2	1	16	0	-8	0	48
                0	1	0	0	1	-15	0	9	0	17
                1	0	0	-2	1	-13	0	7	0	49
                0	-1	0	0	1	-12	0	6	0	19
                2	0	-2	0	0	11	0	0	0	4
                -1	0	2	2	1	-10	0	5	0	50
                1	0	2	2	2	-8	0	3	0	54
                0	-1	2	0	2	-7	0	3	0	53
                0	0	2	2	1	-7	0	3	0	58
                1	1	0	-2	0	-7	0	0	0	51
                0	1	2	0	2	7	0	-3	0	52
                -2	0	0	2	1	-6	0	3	0	20
                0	0	0	2	1	-6	0	3	0	57
                2	0	2	-2	2	6	0	-3	0	56
                1	0	0	2	0	6	0	0	0	55
                1	0	2	-2	1	6	0	-3	0	58
                0	0	0	-2	1	-5	0	3	0	60
                0	-1	2	-2	1	-5	0	3	0	21
                2	0	2	0	1	-5	0	3	0	62
                1	-1	0	0	0	5	0	0	0	61
                1	0	0	-1	0	-4	0	0	0	24
                0	0	0	1	0	-4	0	0	0	65
                0	1	0	-2	0	-4	0	0	0	63
                1	0	-2	0	0	4	0	0	0	64
                2	0	0	-2	1	4	0	-2	0	22
                0	1	2	-2	1	4	0	-2	0	23
                1	1	0	0	0	-3	0	0	0	66
                1	-1	0	-1	0	-3	0	0	0	6
                -1	-1	2	2	2	-3	0	1	0	69
                0	-1	2	2	2	-3	0	1	0	72
                1	-1	2	0	2	-3	0	1	0	68
                3	0	2	0	2	-3	0	1	0	71
                -2	0	2	0	2	-3	0	1	0	5
                1	0	2	0	0	3	0	0	0	67
                -1	0	2	4	2	-2	0	1	0	82
                1	0	0	0	2	-2	0	1	0	76
                -1	0	2	-2	1	-2	0	1	0	74
                0	-2	2	-2	1	-2	0	1	0	7
                -2	0	0	0	1	-2	0	1	0	70
                2	0	0	0	1	2	0	-1	0	75
                3	0	0	0	0	2	0	0	0	77
                1	1	2	0	2	2	0	-1	0	73
                0	0	2	1	2	2	0	-1	0	78
                1	0	0	2	1	-1	0	0	0	91
                1	0	2	2	1	-1	0	1	0	85
                1	1	0	-2	1	-1	0	0	0	102
                0	1	0	2	0	-1	0	0	0	99
                0	1	2	-2	0	-1	0	0	0	30
                0	1	-2	2	0	-1	0	0	0	27
                1	0	-2	2	0	-1	0	0	0	103
                1	0	-2	-2	0	-1	0	0	0	100
                1	0	2	-2	0	-1	0	0	0	94
                1	0	0	-4	0	-1	0	0	0	80
                2	0	0	-4	0	-1	0	0	0	83
                0	0	2	4	2	-1	0	0	0	105
                0	0	2	-1	2	-1	0	0	0	98
                -2	0	2	4	2	-1	0	1	0	86
                2	0	2	2	2	-1	0	0	0	90
                0	-1	2	0	1	-1	0	0	0	101
                0	0	-2	0	1	-1	0	0	0	97
                0	0	4	-2	2	1	0	0	0	92
                0	1	0	0	2	1	0	0	0	28
                1	1	2	-2	2	1	0	-1	0	84
                3	0	2	-2	2	1	0	0	0	93
                -2	0	2	2	2	1	0	-1	0	81
                -1	0	0	0	2	1	0	-1	0	79
                0	0	-2	2	1	1	0	0	0	26
                0	1	2	0	1	1	0	0	0	95
                -1	0	4	0	2	1	0	0	0	87
                2	1	0	-2	0	1	0	0	0	25
                2	0	0	2	0	1	0	0	0	104
                2	0	2	-2	1	1	0	-1	0	89
                2	0	-2	0	1	1	0	0	0	8
                1	-1	0	-2	0	1	0	0	0	88
                -1	0	0	1	1	1	0	0	0	29
                -1	-1	0	2	1	1	0	0	0	96
                0	1	0	1	0	1	0	0	0	106];
        end

        function eopData = getRecentEop(~)
            eopData = [
                2024 10 28 60611  0.217952  0.368278  0.0553462  0.0004272 -0.116789 -0.007716  0.000318  0.000077  37;
                2024 10 29 60612  0.217340  0.367341  0.0548546  0.0005225 -0.116643 -0.007474  0.000316  0.000078  37;
                2024 10 30 60613  0.216738  0.366246  0.0543465  0.0004874 -0.116582 -0.007280  0.000310  0.000076  37;
                2024 10 31 60614  0.216087  0.365657  0.0539091  0.0003834 -0.116573 -0.007282  0.000303  0.000070  37;
                2024 11 01 60615  0.215250  0.364813  0.0535958  0.0002253 -0.116611 -0.007547  0.000298  0.000061  37;
                2024 11 02 60616  0.214668  0.364243  0.0534792  0.0000187 -0.116688 -0.007931  0.000297  0.000053  37;
                2024 11 03 60617  0.213592  0.363739  0.0535507 -0.0001520 -0.116721 -0.008205  0.000299  0.000049  37;
                2024 11 04 60618  0.212708  0.362408  0.0537614 -0.0002789 -0.116601 -0.008276  0.000304  0.000048  37;
                2024 11 05 60619  0.212024  0.361162  0.0540991 -0.0003926 -0.116329 -0.008226  0.000308  0.000048  37;
                2024 11 06 60620  0.211138  0.359897  0.0545311 -0.0004417 -0.116016 -0.008179  0.000313  0.000045  37;
                2024 11 07 60621  0.210020  0.358414  0.0549403 -0.0003609 -0.115747 -0.008168  0.000317  0.000038  37;
                2024 11 08 60622  0.208708  0.356720  0.0552153 -0.0001711 -0.115502 -0.008164  0.000321  0.000029  37;
                2024 11 09 60623  0.207708  0.355248  0.0552937  0.0000617 -0.115238 -0.008148  0.000323  0.000018  37
                2024 11 10 60624  0.206778  0.353838  0.0551146  0.0003232 -0.114978 -0.008096  0.000320  0.000011  37;
                2024 11 11 60625  0.205938  0.352224  0.0546907  0.0005483 -0.114829 -0.007928  0.000310  0.000007  37;
                2024 11 12 60626  0.205036  0.350678  0.0540775  0.0006874 -0.114928 -0.007594  0.000293  0.000007  37;
                2024 11 13 60627  0.204122  0.349304  0.0533965  0.0006785 -0.115282 -0.007236  0.000274  0.000009  37;
                2024 11 14 60628  0.203022  0.348023  0.0528043  0.0005145 -0.115646 -0.007105  0.000259  0.000012  37;
                2024 11 15 60629  0.201871  0.346842  0.0524443  0.0002173 -0.115700 -0.007254  0.000253  0.000016  37;
                2024 11 16 60630  0.200755  0.345545  0.0524216 -0.0001414 -0.115442 -0.007442  0.000257  0.000021  37]';
        end

    end


end