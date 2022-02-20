//Copyright (C) 2021 by Shahar Seifer, Elabum lab, Weizmann Institute of Science
/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace ClusterAlign
{
    public class nogaps
    {

 
        public static double fillgaps(double psi0, double phi0, double[] theta_vec, int S, int contours, int ncenter, int Nx, int Ny, ref int[,] fidx, ref int[,] fidy, ref int[,] fidn, ref double[,] Bfinal, ref double[] Dx_vect, ref double[] Dy_vect, ClusterAlign.Program.tp[,] locations, int[] NFid, int main_IterationNum, ref double Theta_shift, double PreAlignmentTolx, double PreAlignmentToly, ref double xc0, ref double yc0, ref double report_phi, ref double report_psi,string basefilename) //S number of slices
        {
            double pi = Math.PI;
            bool xisrotation = ClusterAlign.Settings4ClusterAlign2.Default.xisRotation; //(Math.Abs(phi0-0.5 * Math.PI)< 0.2 * Math.PI);
            double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = 0;
            double in_a = 0, in_b = 0, in_c = 0, in_d = 0, in_e = 0, in_f = 0, in_g = 0, in_h = 0, in_i = 0;
            double Theta0 = 0;
            double dtheta = theta_vec[2] - theta_vec[1];
            double phi = 0, psi = 0;
            double xc = Nx / 2.0;
            double yc = Ny / 2.0;
             double theta;
            double z0temp1, z0temp2;
            double Dx, Dy;
            double xsim, ysim, zsim;
            double[,] Fx = new double[contours, S];
            double[,] Fy = new double[contours, S];
            double best_angle_psi = psi0;
            double best_angle_phi = phi0;
            double[] zcenter = new double[contours];
            double[] best_zcenter = new double[contours];
            double[] grand_best_zcenter = new double[contours];
            double[] xcenter = new double[contours];
            double[] ycenter = new double[contours];
            double[] previter_xcenter = new double[contours];
            double[] previter_ycenter = new double[contours];
            double[] previter_zcenter = new double[contours];
            double[] grand_Dxbest = new double[ S];
            double[] grand_Dybest = new double[S];
            double[] Dxbest = new double[ S];
            double[] Dybest = new double[S];
            double[] Valt = new double[ S];
            bool[] D_vect_valid = new bool[S];
            bool[] best_D_vect_valid = new bool[ S];
            bool[] grand_best_D_vect_valid = new bool[S];
            Array.Clear(grand_best_D_vect_valid, 0, S);
            double sumerror, error, error2, minerror;
            int count;
            int xc_ns, yc_ns;
            double[] xfactor = new double[S];
            double[] yfactor = new double[S];
            for (int n = 0; n < S; n++)
            {
                if (ClusterAlign.Settings4ClusterAlign2.Default.coswindow)
                {
                    if (xisrotation)
                    {
                        xfactor[n] = 1;
                        yfactor[n] = Math.Cos(theta_vec[n]); //multiplying by these factors means we restore the actual spatial dimensions
                    }
                    else
                    {
                        xfactor[n] = Math.Cos(theta_vec[n]);
                        yfactor[n] = 1;
                    }
                }
                else
                {
                    xfactor[n] = 1;
                    yfactor[n] = 1;
                }
            }

            Console.WriteLine("Fit to rigid body model:");

            //shift to center at xc,yc
            for (int p = 1; p <= contours; p++) 
            {
                for (int ns = 1; ns <= S; ns++)
                {
                    if (ns - 1 == ncenter)
                    {
                        Fx[p - 1, ns - 1] = ((double)fidx[ns - 1, p - 1] - xc0) * xfactor[ns-1];
                        Fy[p - 1, ns - 1] = ((double)fidy[ns - 1, p - 1] - yc0) * yfactor[ns-1];
                    }
                    else
                    {
                        Fx[p - 1, ns - 1] = ((double)fidx[ns - 1, p - 1] - xc) * xfactor[ns - 1];
                        Fy[p - 1, ns - 1] = ((double)fidy[ns - 1, p - 1] - yc) * yfactor[ns - 1];
                    }
                }
            }
            //see below: for nslice: the alignment corrections: Dx=xc-xc0, Dy=yc-yc0,  so fidx,fidy will not change, since it is locations in respect to nonshifted images

            double grand_minerror =0;
            //already defined: double[] Dx_vect = new double[S];
            //already defined: double[] Dy_vect = new double[S];
            double[] Dx_vect_temp = new double[S]; //Dx_vect is updated only if solution found
            double[] Dy_vect_temp = new double[S];

            //using initial values of xcenter, ycenter, zcenter
            for (int p = 1; p <= contours; p++)
            {
                previter_xcenter[p - 1] = Fx[p - 1, ncenter];
                previter_ycenter[p - 1] = Fy[p - 1, ncenter];
                previter_zcenter[p - 1] = 0;
            }



            double[] previter_Dx_vect = new double[S];
            double[] previter_Dy_vect = new double[S];
            //correct for ficticious aspect ratio of images, temporary for this module (then restore back near end)
            for (int ns = 1; ns <= S; ns++)
            {
                Dx_vect[ns - 1] = Dx_vect[ns - 1] * xfactor[ns - 1];
                Dy_vect[ns - 1] = Dy_vect[ns - 1] * yfactor[ns - 1];
            }
            Dx_vect.CopyTo(previter_Dx_vect, 0);//learn from previous rounds
            Dy_vect.CopyTo(previter_Dy_vect, 0);
            

            Array.Clear(zcenter,0, contours); //Start afresh each time, although at the end we publish z0values so matching could improve with iterations
            for (int iteration = 1; iteration <=2; iteration++) //was 4 iterations
            {

                double z_ncenter;
                grand_minerror = 1000000;
                for (phi = (phi0 -15 * pi / 180); phi <= (phi0 + 15 * pi / 180); phi = phi + (1 * pi / 180))//was -15 to 15 step of 3
                {
                    for (psi = (psi0 -10 * pi / 180); psi <= (psi0 + 10* pi / 180); psi = psi + (2 * pi / 180))//tested -20 to 20, step of 4
                    {
                        //for (Theta_shift = -0 * pi / 180; Theta_shift <= 0 * pi / 180; Theta_shift = Theta_shift + 1 * pi / 180)//was step of 2
                        {
                            minerror = 1000000;
                            //double Theta0_best = 0;
                            //for (Theta0 = (-0*0.4 * dtheta + Theta00); Theta0 <= (0*0.4 * dtheta + Theta00); Theta0 = Theta0 + dtheta * 0.2)
                            {

                                if (1==1 ||iteration>1 || main_IterationNum > 0)// 
                                {
                                    // prepare better z0vect
                                    for (int p = 1; p <= contours; p++)
                                    {
                                        int length_vecz0p = 0;
                                        for (int ns = 1; ns <= S; ns++)
                                        {
                                            if (fidn[ns - 1, p - 1] >= 0)
                                            {
                                                length_vecz0p++;
                                            }
                                        }
                                        double z0sum = 0;
                                        int z0count = 0;
                                        double[] vector_z0p = new double[length_vecz0p];
                                        Array.Clear(vector_z0p, 0, length_vecz0p);
                                        for (int ns = 1; ns <= S; ns++)
                                        {
                                            if (fidn[ns-1,p-1]>=0)//((Fx[p - 1, ns - 1] > -xc && Fy[p - 1, ns - 1] > -yc) || ns - 1 == ncenter)
                                            {
                                                theta = theta_vec[ns - 1]  - Theta0;
                                                Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                                                Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                                                z0temp1 = (Fx[p - 1, ns - 1] + previter_Dx_vect[ns - 1]) - a * previter_xcenter[p - 1] - b * previter_ycenter[p - 1];
                                                z0temp2 = (Fy[p - 1, ns - 1] + previter_Dy_vect[ns - 1]) - d * previter_xcenter[p - 1] - e * previter_ycenter[p - 1];
                                                if (Math.Pow(c, 2) + Math.Pow(f, 2) > 0)
                                                {
                                                    vector_z0p[z0count] = (c * z0temp1 + f * z0temp2) / (Math.Pow(c, 2) + Math.Pow(f, 2));
                                                    z0sum = z0sum + (c * z0temp1 + f * z0temp2) / (Math.Pow(c, 2) + Math.Pow(f, 2));
                                                    z0count = z0count + 1;
                                                }
                                            }
                                        }
                                        if (z0count > 0) 
                                        {
                                           zcenter[p - 1] = z0sum / z0count;
                                            if (z0count>8)
                                            {
                                                Array.Sort(vector_z0p);
                                                zcenter[p - 1] = vector_z0p[(int)(length_vecz0p/2)]; //median
                                            }
                                        }
                                    }
                                }

                                theta = Theta0 +  theta_vec[ncenter]; //- to +
                                Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                                Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                                for (int p = 1; p <= contours; p++)
                                {
                                    z_ncenter = i * zcenter[p - 1];
                                    xcenter[p - 1] = in_a * (Fx[p - 1, ncenter]) + in_b * (Fy[p - 1, ncenter]) + in_c * z_ncenter;
                                    ycenter[p - 1] = in_d * (Fx[p - 1, ncenter] )+ in_e * (Fy[p - 1, ncenter]) + in_f * z_ncenter;
                                }
                                // the solved for vairables are in a vector V: Dx[1] Dy[1]..Dx[S] Dy[s]
                                int Vcount,Vsubcount;
                                double Vtempx,Vtempy,Vsumx,Vsumy;
                                Array.Clear(D_vect_valid, 0, S); //default: false

                                for (int ns = 1; ns <= S; ns++)
                                {
                                    if (ns - 1 == ncenter) continue;
                                    int length_vecns = 0;
                                    for (int p = 1; p <= contours; p++)
                                    {
                                        if (fidn[ns - 1, p - 1] >= 0) length_vecns++;
                                    }
                                    double[] vector_dxns = new double[length_vecns];
                                    double[] vector_dyns = new double[length_vecns];
                                    theta = theta_vec[ns - 1] + Theta0;
                                    Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                                    Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                                    Vcount = 0;
                                    Vtempx = 0;
                                    Vtempy = 0;
                                    for (int p = 1; p <= contours; p++)
                                    {
                                        if (fidn[ns - 1, p - 1] >= 0 && fidn[ncenter, p - 1]>=0)
                                        {
                                            vector_dxns[Vcount]= a * xcenter[p - 1] + b * ycenter[p - 1] + c * zcenter[p - 1] - Fx[p - 1, ns - 1];
                                            vector_dyns[Vcount] = d * xcenter[p - 1] + e * ycenter[p - 1] + f * zcenter[p - 1] - Fy[p - 1, ns - 1];
                                            Vtempx = Vtempx + a * xcenter[p - 1] + b * ycenter[p - 1] + c * zcenter[p - 1] - Fx[p - 1, ns - 1];
                                            Vtempy = Vtempy + d * xcenter[p - 1] + e * ycenter[p - 1] + f * zcenter[p - 1] - Fy[p - 1, ns - 1];
                                            Vcount = Vcount + 1;

                                        }
                                    }
                                    if (Vcount > 0)
                                    {
                                        Dx_vect_temp[ns - 1] = Vtempx / Vcount;
                                        Dy_vect_temp[ns - 1] = Vtempy / Vcount;
                                        D_vect_valid[ns - 1] = true;
                                         if (Vcount>=9)
                                        {
                                            Array.Sort(vector_dxns);
                                            Array.Sort(vector_dyns);
                                            Vsumx = 0;
                                            Vsumy = 0;
                                            Vsubcount = 0;
                                            for (int t=(int)(0.4*length_vecns);t< length_vecns-(int)(0.4 * length_vecns); t++)
                                            {
                                                Vsumx+= vector_dxns[t];
                                                Vsumy += vector_dyns[t];
                                                Vsubcount++;
                                            }
                                            Dx_vect_temp[ns - 1] = Vsumx / Vsubcount;
                                            Dy_vect_temp[ns - 1] = Vsumy / Vsubcount;
                                        }
                                        if (Math.Abs(Dx_vect_temp[ns - 1]) > PreAlignmentTolx)
                                        {
                                            Dx_vect_temp[ns - 1] = Math.Sign(Dx_vect_temp[ns - 1]) * PreAlignmentTolx;
                                        }
                                        if (Math.Abs(Dy_vect_temp[ns - 1]) > PreAlignmentToly)
                                        {
                                            Dy_vect_temp[ns - 1] = Math.Sign(Dy_vect_temp[ns - 1]) * PreAlignmentToly;
                                        }
                                    }
                                    else
                                    {
                                        D_vect_valid[ns - 1] = false;
                                        Dx_vect_temp[ns - 1] = previter_Dx_vect[ns - 1];
                                        Dy_vect_temp[ns - 1] = previter_Dy_vect[ns - 1];
                                    }
                                }
                                sumerror = 0;
                                count = 0;
                                for (int p = 1; p <= contours; p++)
                                {
                                    for (int ns = 1; ns <= S; ns++)
                                    {
                                        theta = theta_vec[ns - 1] + Theta_shift - Theta0;
                                        Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                                        Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                                        if (fidn[ns - 1, p - 1] >= 0 && fidn[ncenter, p - 1] >= 0)
                                        {
                                            xsim = a * xcenter[p - 1] + b * ycenter[p - 1] + c * zcenter[p - 1]; // z0vect(p) replaced with initial zcenter
                                            ysim = d * xcenter[p - 1] + e * ycenter[p - 1] + f * zcenter[p - 1];
                                            zsim = g * xcenter[p - 1] + h * ycenter[p - 1] + i * zcenter[p - 1];
                                            error2 = Math.Pow((Fx[p - 1, ns - 1] + Dx_vect_temp[ns - 1] - xsim), 2) + Math.Pow((Fy[p - 1, ns - 1] + Dy_vect_temp[ns - 1] - ysim), 2);
                                            sumerror = sumerror + error2;
                                            count = count + 1;
                                        }
                                    }
                                }
                                error = Math.Sqrt(sumerror / count);
                                if (error < minerror)
                                {
                                    minerror = error;
                                    zcenter.CopyTo(best_zcenter, 0); //second argument is the starting index to copy (length: contour)
                                    Dx_vect_temp.CopyTo(Dxbest, 0); //length S
                                    Dy_vect_temp.CopyTo(Dybest, 0); //length S
                                    D_vect_valid.CopyTo(best_D_vect_valid, 0);
                                }
                            }
                            if (minerror < grand_minerror)
                            {
                                grand_minerror = minerror;
                                best_angle_psi = psi;
                                best_angle_phi = phi;
                                Dxbest.CopyTo(grand_Dxbest, 0);
                                Dybest.CopyTo(grand_Dybest, 0);
                                best_D_vect_valid.CopyTo(grand_best_D_vect_valid, 0);
                                 best_zcenter.CopyTo(grand_best_zcenter, 0);
                            }
                        } //for Theta_shift
                    } //for psi
                } //for phi


                psi = best_angle_psi;
                phi = best_angle_phi;
                for (int tempn = 1; tempn <= S; tempn = tempn + 1)
                {
                    if (grand_best_D_vect_valid[tempn - 1] && Math.Abs(grand_Dxbest[tempn - 1]) <= PreAlignmentTolx && Math.Abs(grand_Dybest[tempn - 1]) <= PreAlignmentToly)
                    {
                        previter_Dx_vect[tempn - 1] = grand_Dxbest[tempn - 1];
                        previter_Dy_vect[tempn - 1] = grand_Dybest[tempn - 1];
                    }
                 }

                grand_best_zcenter.CopyTo(zcenter,0);

                theta = Theta0 + theta_vec[ncenter];//- to +
                Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                for (int p = 1; p <= contours; p++)
                {
                    z_ncenter = i * zcenter[p - 1];
                    previter_xcenter[p - 1] = in_a * Fx[p - 1, ncenter] + in_b * Fy[p - 1, ncenter] + in_c * z_ncenter;
                    previter_ycenter[p - 1] = in_d * Fx[p - 1, ncenter] + in_e * Fy[p - 1, ncenter] + in_f * z_ncenter;
                    previter_zcenter[p - 1] = zcenter[p - 1];
                }

 



            } //for iteration

            Console.WriteLine("phi={0:0.0}, psi={1:0.0}", phi * 180 / pi,psi * 180 / pi);
            report_phi = phi;
            report_psi = psi;


            previter_xcenter.CopyTo(xcenter, 0);
            previter_ycenter.CopyTo(ycenter, 0);
            previter_zcenter.CopyTo(zcenter, 0);

            //Build or correct Dx_vect, Dy_vect (planar movement compensations): in case that the values are still not valid use data of all fiducials at center extrapoalted by model to other slices compared to actual points
            //if values in previous iteration are valid then use those without additonal calculations.
            //Fill_DxDy(previter_Dx_vect, previter_Dy_vect, grand_best_D_vect_valid, NFid, locations, Nx, Ny, xcenter, ycenter, zcenter, S, contours, theta_vec, Theta_shift, Theta0, phi, psi, Dx_vect, Dy_vect);
            previter_Dx_vect.CopyTo(Dx_vect, 0);//store results in permanent vectors
            previter_Dy_vect.CopyTo(Dy_vect, 0);

            //write_testfile(ref z0vect, "D:\\results\\z0vect.txt");
            write_testfile(ref xcenter, basefilename + ".xcenter.txt");
            write_testfile(ref ycenter, basefilename + ".ycenter.txt");
            write_testfile(ref zcenter, basefilename + ".zcenter.txt");
            write_testfile(ref Dx_vect, basefilename + ".Dx.txt");
            write_testfile(ref Dy_vect, basefilename + ".Dy.txt");

            //correct back to aspect ratio of images
            for (int ns = 1; ns <= S; ns++)
            {
                Dx_vect[ns - 1] = Dx_vect[ns - 1]/xfactor[ns-1];
                Dy_vect[ns - 1] = Dy_vect[ns - 1]/yfactor[ns-1];
            }
            //Fill Bfinal table
            int cont, slice;
            for (int n = 0; n <= S * contours - 1; n++)
            {
                cont = (int)(n / S);
                slice = n % S;
                if (fidn[slice, cont] >= 0 || slice == ncenter)//((Fx[cont + 1 - 1, slice + 1 - 1] > -xc && Fy[cont + 1 - 1, slice + 1 - 1] > -yc) || slice==ncenter)
                {
                     xc_ns = (int)(slice == ncenter ? xc0 : xc);
                    yc_ns = (int)(slice == ncenter ? yc0 : yc);
                    Bfinal[n + 1 - 1, 1 - 1] = 0;
                    Bfinal[n + 1 - 1, 2 - 1] = cont;
                    Bfinal[n + 1 - 1, 3 - 1] = Math.Round(Fx[cont + 1 - 1, slice + 1 - 1] / xfactor[slice] + xc_ns);
                    Bfinal[n + 1 - 1, 4 - 1] = Math.Round(Fy[cont + 1 - 1, slice + 1 - 1] / yfactor[slice] + yc_ns);
                    Bfinal[n + 1 - 1, 5 - 1] = slice;
                }
            }
            double[] fit_err = new double[S];
            for (int ns = 1; ns <= S; ns++)
            {
                theta = theta_vec[ns - 1] + Theta_shift - Theta0;
                Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                Dx = Dx_vect[ns - 1];
                Dy = Dy_vect[ns - 1];
                sumerror = 0;
                count = 0;
                for (int p = 1; p <= contours; p++)
                {
                    xsim = a * xcenter[p - 1] + b * ycenter[p - 1] + c * zcenter[p - 1];
                    ysim = d * xcenter[p - 1] + e * ycenter[p - 1] + f * zcenter[p - 1];
                    zsim = g * xcenter[p - 1] + h * ycenter[p - 1] + i * zcenter[p - 1];
                    if (fidn[ns - 1, p - 1]<0 || fidn[ncenter, p - 1] < 0)
                    {
                        // implement simulation to fill gaps
                        Bfinal[(p - 1) * S + ns - 1, 1 - 1] = 0;
                        Bfinal[(p - 1) * S + ns - 1, 2 - 1] = p - 1; //contour (fiducial) number
                        Bfinal[(p - 1) * S + ns - 1, 3 - 1] = Math.Round(xsim / xfactor[ns-1] - Dx + xc);
                        Bfinal[(p - 1) * S + ns - 1, 4 - 1] = Math.Round(ysim / yfactor[ns-1] - Dy + yc);
                        Bfinal[(p - 1) * S + ns - 1, 5 - 1] = ns - 1; //slice number
                        //fill back the fiducial tables in program.c
                        fidx[ns - 1, p - 1] = (int)Math.Round(xsim / xfactor[ns - 1] - Dx + xc);
                        fidy[ns - 1, p - 1] = (int)Math.Round(ysim / yfactor[ns - 1] - Dy + yc);
                        fidn[ns - 1, p - 1] = -2; //still mark it as non-authentic point. Otherwise would be =p-1;
                    }
                    else //if points of FX/Fy do exit, then calculate error
                    {
                        error2 = Math.Pow((Fx[p - 1, ns - 1] + Dx_vect[ns - 1] - xsim), 2) + Math.Pow((Fy[p - 1, ns - 1] + Dy_vect[ns - 1] - ysim), 2);
                        sumerror = sumerror + error2;
                        count = count + 1;
                    }
                }
                if (count > 0)
                { fit_err[ns - 1] = Math.Sqrt(sumerror / count); }
                else
                { fit_err[ns - 1] = -1; }
            }
            write_testfile(ref fit_err, basefilename + ".fit_err_by_slice.txt");

            return grand_minerror;

        }



        public static void Amat(double theta, double phi, double psi, ref double a, ref double b, ref double c, ref double d, ref double e, ref double f, ref double g, ref double h, ref double i) //provides all elements of Amat
        {
            a = Math.Cos(theta) + (1 - Math.Cos(theta)) * Math.Pow(Math.Cos(psi), 2) * Math.Pow(Math.Sin(phi), 2);
            b = 0 + (1 - Math.Cos(theta)) * Math.Pow(Math.Cos(psi), 2) * Math.Cos(phi) * Math.Sin(phi) - Math.Sin(theta) * Math.Sin(psi);
            c = 0 + (1 - Math.Cos(theta)) * Math.Cos(psi) * Math.Sin(psi) * Math.Sin(phi) + Math.Sin(theta) * Math.Cos(psi) * Math.Cos(phi);
            d = 0 + (1 - Math.Cos(theta)) * Math.Pow(Math.Cos(psi), 2) * Math.Cos(phi) * Math.Sin(phi) + Math.Sin(theta) * Math.Sin(psi);
            e = Math.Cos(theta) + (1 - Math.Cos(theta)) * Math.Pow(Math.Cos(psi), 2) * Math.Pow(Math.Cos(phi), 2);
            f = 0 + (1 - Math.Cos(theta)) * Math.Cos(psi) * Math.Sin(psi) * Math.Cos(phi) - Math.Sin(theta) * Math.Cos(psi) * Math.Sin(phi);
            g = 0 + (1 - Math.Cos(theta)) * Math.Cos(psi) * Math.Sin(psi) * Math.Sin(phi) - Math.Sin(theta) * Math.Cos(psi) * Math.Cos(phi);
            h = 0 + (1 - Math.Cos(theta)) * Math.Cos(psi) * Math.Sin(psi) * Math.Cos(phi) + Math.Sin(theta) * Math.Cos(psi) * Math.Sin(phi);
            i = Math.Cos(theta) + (1 - Math.Cos(theta)) * Math.Pow(Math.Sin(psi), 2);
            double thresh = 1e-6;
            if (Math.Abs(a) < thresh) a = 0;
            if (Math.Abs(b) < thresh) b = 0;
            if (Math.Abs(c) < thresh) c = 0;
            if (Math.Abs(d) < thresh) d = 0;
            if (Math.Abs(e) < thresh) e = 0;
            if (Math.Abs(f) < thresh) f = 0;
            if (Math.Abs(g) < thresh) g = 0;
            if (Math.Abs(h) < thresh) h = 0;
            if (Math.Abs(i) < thresh) i = 0;
        }

        public static void Inv(double a, double b, double c, double d, double e, double f, double g, double h, double i, ref double in_a, ref double in_b, ref double in_c, ref double in_d, ref double in_e, ref double in_f, ref double in_g, ref double in_h, ref double in_i)
        {
            double G = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
            in_a = (e * i - f * h) / G;
            in_b = (c * h - b * i) / G;
            in_c = (b * f - c * e) / G;
            in_d = (f * g - d * i) / G;
            in_e = (a * i - c * g) / G;
            in_f = (c * d - a * f) / G;
            in_g = (d * h - e * g) / G;
            in_h = (b * g - a * h) / G;
            in_i = (a * e - b * d) / G;
        }

        public static void write_testfile(ref double[] field, string FileName)
        {
            var sr = new StreamWriter(FileName, true, System.Text.Encoding.UTF8);//true: append
            //sr.NewLine = "\n";
            string gap="";
            for (int i = 0; i < field.Length; i++)
            {
                if (i < field.Length - 1)
                {gap = "\t"; }
                else
                { gap = "\n"; }
                sr.Write(field[i].ToString() + gap);
            }
            sr.Close();
        }


        public static double[] estimate_z(bool xisrotation, int ncenter, int Nx, int Ny, ClusterAlign.Program.tp[,] locations, int[] NFid, double PreAlignmentTolx, double PreAlignmentToly, ref double[] tiltangles)
        //OLD code, not used, not effective
        {
            double maxz = 1000;
            double fit_tol = 5;
            int pstep = 3, mstep = 3;
            double xc = Nx / 2.0;
            double yc = Ny / 2.0;
            int N = NFid[ncenter];
            double[] z_values = new double[N];
            Array.Clear(z_values, 0, N);
            //find zvalues for fiducials in ncenter slice based on correlation with locations at ncenter+1 and ncenter-1
            //step1: tentative image shifts
            int shiftx_p = 0, shiftx_m = 0, shifty_p = 0, shifty_m = 0;
            detect_DxDy_mp(NFid, locations, Nx, Ny, ncenter, mstep, pstep, N, PreAlignmentTolx, PreAlignmentToly, ref shiftx_p, ref shiftx_m, ref shifty_p, ref shifty_m);
            //step2: choose the most likely track (near target + fitting to consistent model)
            double costhm = Math.Cos(tiltangles[ncenter - mstep]);
            double costhc = Math.Cos(tiltangles[ncenter]);
            double costhp = Math.Cos(tiltangles[ncenter + pstep]);
            double sinthm = Math.Sin(tiltangles[ncenter - mstep]);
            double sinthc = Math.Sin(tiltangles[ncenter]);
            double sinthp = Math.Sin(tiltangles[ncenter + pstep]);
            double rangem = Math.Abs(maxz * (sinthm - sinthc));
            double rangep = Math.Abs(maxz * (sinthp - sinthc));
            double x, y, xm, ym, xp, yp, targetxm, targetxp, targetym, targetyp;
            double delta, line_err, min_distance;
            for (int nc = 0; nc < NFid[ncenter]; nc++)
            {
                x = locations[ncenter, nc].col - xc;
                y = locations[ncenter, nc].row - yc;
                if (xisrotation)
                {
                    targetxm = x + shiftx_m; //still no clue about z, it's the best location to expect the fiducial in slice ncenter-1
                    targetym = y * costhm / costhc + shifty_m;
                }
                else
                {
                    targetxm = x * costhm / costhc + shiftx_m; //still no clue about z, it's the best location to expect the fiducial in slice ncenter-1
                    targetym = y + shifty_m;
                }
                min_distance = rangem;
                for (int nm = 0; nm < NFid[ncenter - mstep]; nm++)
                {
                    xm = locations[ncenter - mstep, nm].col - xc;
                    ym = locations[ncenter - mstep, nm].row - yc;
                    if (xisrotation)
                    {
                        delta = (ym - targetym);
                        line_err = Math.Abs(xm - targetxm);
                    }
                    else
                    {
                        delta = (xm - targetxm);
                        line_err = Math.Abs(ym - targetym);
                    }
                    if (Math.Abs(delta) < min_distance && line_err <= fit_tol)
                    {
                        if (xisrotation)
                        {
                            targetxp = x + shiftx_p;
                            targetyp = y * costhp / costhc + delta * (sinthp - sinthc) / (sinthm - sinthc) + shifty_p; //delta expected to be z*(sinthm-sinthc)
                        }
                        else
                        {
                            targetxp = x * costhp / costhc + delta * (sinthp - sinthc) / (sinthm - sinthc) + shiftx_p; //delta expected to be z*(sinthm-sinthc)
                            targetyp = y + shifty_p;
                        }

                        for (int np = 0; np < NFid[ncenter + pstep]; np++)
                        {
                            xp = locations[ncenter + pstep, np].col - xc;
                            yp = locations[ncenter + pstep, np].row - yc;
                            if (Math.Abs(xp - targetxp) <= fit_tol && Math.Abs(yp - targetyp) <= fit_tol)
                            {
                                min_distance = Math.Abs(delta);
                                z_values[nc] = delta / (sinthm - sinthc);
                            }

                        }
                    }
                }
            }
            //write_testfile(ref z_values, "D:\\results\\z_initial.txt");
            return z_values;
        }

        public static void detect_DxDy_mp(int[] NFid, ClusterAlign.Program.tp[,] locations, int Nx, int Ny, int ncenter, int mstep, int pstep, int contours, double PreAlignmentTolx, double PreAlignmentToly, ref int shiftx_p, ref int shiftx_m, ref int shifty_p, ref int shifty_m)
        //OLD code, not used and less effective
        {
            //Compare ncenter image to ncenter-1 and ncenter+1 and derive image shifts 
            int coarse_factor = 8;
            int bin = 4;
            int range_factor = coarse_factor;
            double xc = Nx / 2.0;
            double yc = Ny / 2.0;
            int nx_coarse = (int)((Nx / bin) / coarse_factor) + 1;
            int ny_coarse = (int)((Ny / bin) / coarse_factor) + 1;
            bool[,] mesh1 = new bool[nx_coarse, ny_coarse];
            bool[,] mesh2 = new bool[nx_coarse, ny_coarse];
            bool[,] mesh1_full = new bool[nx_coarse * coarse_factor, ny_coarse * coarse_factor];
            bool[,] mesh2_full = new bool[nx_coarse * coarse_factor, ny_coarse * coarse_factor];
            int count_mesh;
            int count_mesh_max;
            int cost, cost_max, cost_fine, cost_fine_max;
            int shift_Dx_fine = 0, shift_Dy_fine = 0, shift_Dx_fine_best = 0, shift_Dy_fine_best = 0;
            int sh_limit = (int)(nx_coarse * 0.3);
            Int64 sum_shift_Dx;
            Int64 sum_shift_Dy;
            int xpnt, ypnt;
            int xpnt_coarse, ypnt_coarse, xpnt_coarse1, ypnt_coarse1;
            int shift_Dx_best = 0, shift_Dy_best = 0;
            for (int ns = (ncenter- mstep) +1; ns <= (ncenter + pstep) + 1; ns=ns+mstep+pstep)
            {
                {
                    Console.WriteLine("Slice {0} alignment is compromised due to low number of tracked fiducials.", ns);
                    //Overwrite values of Dx_vect, Dy_vect in slices ncenter-1 ->m,  and ncenter+1 ->p 
                    //prepare mesh1 with fiducial points in slice ncenter
                    Array.Clear(mesh1, 0, nx_coarse * ny_coarse);
                    Array.Clear(mesh1_full, 0, nx_coarse * coarse_factor * ny_coarse * coarse_factor);
                    int count_points = 0;
                    for (int p = 1; p <= contours; p++)
                    {
                        xpnt = locations[ncenter, p - 1].col;
                        ypnt = locations[ncenter, p - 1].row;
                        if (xpnt >= 0 && xpnt < Nx && ypnt >= 0 && ypnt < Ny)
                        {
                            xpnt_coarse = Math.Max(Math.Min((int)((xpnt / bin) / coarse_factor), nx_coarse - 1), 0);
                            ypnt_coarse = Math.Max(Math.Min((int)((ypnt / bin) / coarse_factor), ny_coarse - 1), 0);
                            mesh1[xpnt_coarse, ypnt_coarse] = true;
                            count_points++;
                            mesh1_full[(int)(xpnt / bin), (int)(ypnt / bin)] = true;
                        }
                    }
                    //prepare mesh2 with information in locations (all optically reconginzed fiducials)
                    Array.Clear(mesh2, 0, nx_coarse * ny_coarse);
                    Array.Clear(mesh2_full, 0, nx_coarse * coarse_factor * ny_coarse * coarse_factor);
                    for (int pnt = 1; pnt <= NFid[ns - 1]; pnt++)
                    {
                        xpnt = locations[ns - 1, pnt - 1].col;
                        ypnt = locations[ns - 1, pnt - 1].row;
                        if (xpnt < 0 || ypnt < 0 || xpnt > Nx - 1 || ypnt > Ny - 1) continue;
                        mesh2_full[(int)(xpnt / bin), (int)(ypnt / bin)] = true;
                        xpnt_coarse = (int)((xpnt / bin) / coarse_factor);
                        ypnt_coarse = (int)((ypnt / bin) / coarse_factor);
                        mesh2[xpnt_coarse, ypnt_coarse] = true;
                    }
                    count_mesh_max = 0;
                    cost_max = 0;
                    for (int shift_Dx = -sh_limit; shift_Dx <= sh_limit; shift_Dx++)
                    {
                        for (int shift_Dy = -sh_limit; shift_Dy <= sh_limit; shift_Dy++)
                        {
                            cost = 0;
                            count_mesh = 0;
                            sum_shift_Dx = 0;
                            sum_shift_Dy = 0;
                            int xpnt_coarse_min = Math.Max(0, -shift_Dx);
                            int ypnt_coarse_min = Math.Max(0, -shift_Dy);
                            int xpnt_coarse_max = Math.Min(nx_coarse, nx_coarse - shift_Dx);
                            int ypnt_coarse_max = Math.Min(ny_coarse, ny_coarse - shift_Dy);

                            for (xpnt_coarse = xpnt_coarse_min; xpnt_coarse < xpnt_coarse_max; xpnt_coarse++)
                            {
                                for (ypnt_coarse = ypnt_coarse_min; ypnt_coarse < ypnt_coarse_max; ypnt_coarse++)
                                {
                                    xpnt_coarse1 = xpnt_coarse + shift_Dx;
                                    ypnt_coarse1 = ypnt_coarse + shift_Dy;
                                    if (mesh1[xpnt_coarse1, ypnt_coarse1] && mesh2[xpnt_coarse, ypnt_coarse])
                                    {
                                        cost_fine_max = 0;
                                        shift_Dx_fine_best = 0;
                                        shift_Dy_fine_best = 0;
                                        for (shift_Dx_fine = -range_factor + 1; shift_Dx_fine < range_factor; shift_Dx_fine++)
                                        {
                                            for (shift_Dy_fine = -range_factor + 1; shift_Dy_fine < range_factor; shift_Dy_fine++)
                                            {
                                                int min_tempx = Math.Max(0, -shift_Dx_fine - xpnt_coarse1 * coarse_factor);
                                                int max_tempx = Math.Min(coarse_factor, nx_coarse * coarse_factor - xpnt_coarse1 * coarse_factor - shift_Dx_fine);
                                                int min_tempy = Math.Max(0, -shift_Dy_fine - ypnt_coarse1 * coarse_factor);
                                                int max_tempy = Math.Min(coarse_factor, ny_coarse * coarse_factor - ypnt_coarse1 * coarse_factor - shift_Dy_fine);
                                                int ref1x = shift_Dx_fine + xpnt_coarse1 * coarse_factor;
                                                int ref1y = shift_Dy_fine + ypnt_coarse1 * coarse_factor;
                                                int ref2x = xpnt_coarse * coarse_factor;
                                                int ref2y = ypnt_coarse * coarse_factor;
                                                cost_fine = 0;
                                                for (int tempx = min_tempx; tempx < max_tempx; tempx++)
                                                {
                                                    for (int tempy = min_tempy; tempy < max_tempy; tempy++)
                                                    {
                                                        //correlation is positive if matching shift distance is up to range_factor pixels
                                                        if (mesh1_full[ref1x + tempx, ref1y + tempy] && mesh2_full[ref2x + tempx, ref2y + tempy])
                                                        {
                                                            cost_fine++;
                                                        }
                                                    }
                                                }
                                                if (cost_fine > cost_fine_max)
                                                {
                                                    cost_fine_max = cost_fine;
                                                    shift_Dx_fine_best = shift_Dx_fine;
                                                    shift_Dy_fine_best = shift_Dy_fine;
                                                }
                                                //if (flag) break;
                                            }
                                            //if (flag) break;
                                        }
                                        if (cost_fine_max > 0)
                                        {
                                            count_mesh++;
                                            cost = cost + cost_fine_max;
                                            sum_shift_Dx = sum_shift_Dx + shift_Dx * coarse_factor + shift_Dx_fine_best;
                                            sum_shift_Dy = sum_shift_Dy + shift_Dy * coarse_factor + shift_Dy_fine_best;
                                        }
                                    }
                                }
                            }
                            if (cost > cost_max)
                            {
                                cost_max = cost;
                                count_mesh_max = count_mesh;
                                shift_Dx_best = bin * (int)(sum_shift_Dx / count_mesh);
                                shift_Dy_best = bin * (int)(sum_shift_Dy / count_mesh);
                            }
                        }
                    }
                    if (count_mesh_max > 0.5 * count_points && Math.Abs(shift_Dx_best) <= 5*PreAlignmentTolx && Math.Abs(shift_Dy_best) <= 5*PreAlignmentToly)
                    {
                        if (ns-1==ncenter-mstep)
                        {
                            shiftx_m = shift_Dx_best;
                            shifty_m = shift_Dy_best;
                        }
                        else
                        {
                            shiftx_p = shift_Dx_best;
                            shifty_p = shift_Dy_best;
                        }
                    }
                }
            }

        }



        public static void Fill_DxDy(bool[] grand_best_D_vect_valid, int[] NFid, ClusterAlign.Program.tp[,] locations, int Nx, int Ny, double[] xcenter, double[] ycenter, double[] zcenter, int S, int contours, double[] theta_vec, double Theta_shift, double Theta0, double phi, double psi, double[] Dx_vect, double[] Dy_vect, double PreAlignmentTolx, double PreAlignmentToly)
        //OLD code, not used and less effective
        {
            int coarse_factor = 8;
            int bin = 4;
            int range_factor = coarse_factor;
            double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = 0;
            double in_a = 0, in_b = 0, in_c = 0, in_d = 0, in_e = 0, in_f = 0, in_g = 0, in_h = 0, in_i = 0;
            double xc = Nx / 2.0;
            double yc = Ny / 2.0;
            int nx_coarse = (int)((Nx / bin) / coarse_factor) + 1;
            int ny_coarse = (int)((Ny / bin) / coarse_factor) + 1;
            bool[,] mesh1 = new bool[nx_coarse, ny_coarse];
            bool[,] mesh2 = new bool[nx_coarse, ny_coarse];
            bool[,] mesh1_full = new bool[nx_coarse * coarse_factor, ny_coarse * coarse_factor];
            bool[,] mesh2_full = new bool[nx_coarse * coarse_factor, ny_coarse * coarse_factor];
            double theta, Dx, Dy, xsim, ysim;
            int count_mesh;
            int count_mesh_max;
            int cost, cost_max, cost_fine, cost_fine_max;
            int shift_Dx_fine = 0, shift_Dy_fine = 0, shift_Dx_fine_best = 0, shift_Dy_fine_best = 0;
            int sh_limit = (int)(nx_coarse * 0.3);
            Int64 sum_shift_Dx;
            Int64 sum_shift_Dy;
            int xpnt, ypnt;
            int xpnt_coarse, ypnt_coarse, xpnt_coarse1, ypnt_coarse1;
            int shift_Dx_best = 0, shift_Dy_best = 0;
            for (int ns = 1; ns <= S; ns++)
            {
                if (grand_best_D_vect_valid[(ns - 1)*2] == false || grand_best_D_vect_valid[(ns - 1) * 2+1] == false)
                {
                    Console.WriteLine("Slice {0} alignment is compromised due to low number of tracked fiducials.", ns);
                    //Overwrite missing/wrong values of Dx_vect, Dy_vect according to all fiducial fitting to simulated values 
                    //prepare mesh1 with simulated fidx,fidy points completly based on model to avoid confusion with multiple Dx/Dy values
                    Array.Clear(mesh1, 0, nx_coarse * ny_coarse);
                    Array.Clear(mesh1_full, 0, nx_coarse * coarse_factor * ny_coarse * coarse_factor);
                    int count_points = 0;
                    theta = theta_vec[ns - 1]+ Theta_shift - Theta0;
                    Amat(theta, phi, psi, ref a, ref b, ref c, ref d, ref e, ref f, ref g, ref h, ref i);
                    Inv(a, b, c, d, e, f, g, h, i, ref in_a, ref in_b, ref in_c, ref in_d, ref in_e, ref in_f, ref in_g, ref in_h, ref in_i);
                    Dx = Dx_vect[ns - 1];
                    Dy = Dy_vect[ns - 1];
                    for (int p = 1; p <= contours; p++)
                    {
                        xsim = a * xcenter[p - 1] + b * ycenter[p - 1] + c * zcenter[p - 1];
                        ysim = d * xcenter[p - 1] + e * ycenter[p - 1] + f * zcenter[p - 1];
                        xpnt = (int)Math.Round(xsim - Dx + xc);
                        ypnt = (int)Math.Round(ysim - Dy + yc);
                        if (xpnt >= 0 && xpnt < Nx && ypnt >= 0 && ypnt < Ny)
                        {
                            xpnt_coarse = Math.Max(Math.Min((int)((xpnt / bin) / coarse_factor), nx_coarse - 1), 0);
                            ypnt_coarse = Math.Max(Math.Min((int)((ypnt / bin) / coarse_factor), ny_coarse - 1), 0);
                            mesh1[xpnt_coarse, ypnt_coarse] = true;
                            count_points++;
                            mesh1_full[(int)(xpnt / bin), (int)(ypnt / bin)] = true;
                        }
                    }
                    //prepare mesh2 with information in locations (all optically reconginzed fiducials)
                    Array.Clear(mesh2, 0, nx_coarse * ny_coarse);
                    Array.Clear(mesh2_full, 0, nx_coarse * coarse_factor * ny_coarse * coarse_factor);
                    for (int pnt = 1; pnt <= NFid[ns - 1]; pnt++)
                    {
                        xpnt = locations[ns - 1, pnt - 1].col;
                        ypnt = locations[ns - 1, pnt - 1].row;
                        if (xpnt < 0 || ypnt < 0 || xpnt > Nx - 1 || ypnt > Ny - 1) continue;
                        mesh2_full[(int)(xpnt / bin), (int)(ypnt / bin)] = true;
                        xpnt_coarse = (int)((xpnt / bin) / coarse_factor);
                        ypnt_coarse = (int)((ypnt / bin) / coarse_factor);
                        mesh2[xpnt_coarse, ypnt_coarse] = true;
                    }
                    count_mesh_max = 0;
                    cost_max = 0;
                    for (int shift_Dx = -sh_limit; shift_Dx <= sh_limit; shift_Dx++)
                    {
                        for (int shift_Dy = -sh_limit; shift_Dy <= sh_limit; shift_Dy++)
                        {
                            cost = 0;
                            count_mesh = 0;
                            sum_shift_Dx = 0;
                            sum_shift_Dy = 0;
                            int xpnt_coarse_min = Math.Max(0, -shift_Dx);
                            int ypnt_coarse_min = Math.Max(0, -shift_Dy);
                            int xpnt_coarse_max = Math.Min(nx_coarse, nx_coarse - shift_Dx);
                            int ypnt_coarse_max = Math.Min(ny_coarse, ny_coarse - shift_Dy);

                            for (xpnt_coarse = xpnt_coarse_min; xpnt_coarse < xpnt_coarse_max; xpnt_coarse++)
                            {
                                for (ypnt_coarse = ypnt_coarse_min; ypnt_coarse < ypnt_coarse_max; ypnt_coarse++)
                                {
                                    xpnt_coarse1 = xpnt_coarse + shift_Dx;
                                    ypnt_coarse1 = ypnt_coarse + shift_Dy;
                                    if (mesh1[xpnt_coarse1, ypnt_coarse1] && mesh2[xpnt_coarse, ypnt_coarse])
                                    {
                                        cost_fine_max = 0;
                                        shift_Dx_fine_best = 0;
                                        shift_Dy_fine_best = 0;
                                        for (shift_Dx_fine = -range_factor + 1; shift_Dx_fine < range_factor; shift_Dx_fine++)
                                        {
                                            for (shift_Dy_fine = -range_factor + 1; shift_Dy_fine < range_factor; shift_Dy_fine++)
                                            {
                                                int min_tempx = Math.Max(0, -shift_Dx_fine - xpnt_coarse1 * coarse_factor);
                                                int max_tempx = Math.Min(coarse_factor, nx_coarse * coarse_factor - xpnt_coarse1 * coarse_factor - shift_Dx_fine);
                                                int min_tempy = Math.Max(0, -shift_Dy_fine - ypnt_coarse1 * coarse_factor);
                                                int max_tempy = Math.Min(coarse_factor, ny_coarse * coarse_factor - ypnt_coarse1 * coarse_factor - shift_Dy_fine);
                                                int ref1x = shift_Dx_fine + xpnt_coarse1 * coarse_factor;
                                                int ref1y = shift_Dy_fine + ypnt_coarse1 * coarse_factor;
                                                int ref2x = xpnt_coarse * coarse_factor;
                                                int ref2y = ypnt_coarse * coarse_factor;
                                                cost_fine = 0;
                                                for (int tempx = min_tempx; tempx < max_tempx; tempx++)
                                                {
                                                    for (int tempy = min_tempy; tempy < max_tempy; tempy++)
                                                    {
                                                        //correlation is positive if matching shift distance is up to range_factor pixels
                                                        if (mesh1_full[ref1x + tempx, ref1y + tempy] && mesh2_full[ref2x + tempx, ref2y + tempy])
                                                        {
                                                            cost_fine++;
                                                        }
                                                    }
                                                }
                                                if (cost_fine > cost_fine_max)
                                                {
                                                    cost_fine_max = cost_fine;
                                                    shift_Dx_fine_best = shift_Dx_fine;
                                                    shift_Dy_fine_best = shift_Dy_fine;
                                                }
                                                //if (flag) break;
                                            }
                                            //if (flag) break;
                                        }
                                        if (cost_fine_max > 0)
                                        {
                                            count_mesh++;
                                            cost = cost + cost_fine_max;
                                            sum_shift_Dx = sum_shift_Dx + shift_Dx * coarse_factor + shift_Dx_fine_best;
                                            sum_shift_Dy = sum_shift_Dy + shift_Dy * coarse_factor + shift_Dy_fine_best;
                                        }
                                    }
                                }
                            }
                            if (cost > cost_max)
                            {
                                cost_max = cost;
                                count_mesh_max = count_mesh;
                                shift_Dx_best = bin * (int)(sum_shift_Dx / count_mesh);
                                shift_Dy_best = bin * (int)(sum_shift_Dy / count_mesh);
                            }
                        }
                    }
                    if (count_mesh_max > 0.75 * count_points &&  Math.Abs(Dx_vect[ns - 1] + shift_Dx_best )<= PreAlignmentTolx && Math.Abs(Dy_vect[ns - 1] + shift_Dy_best) <= PreAlignmentToly)
                    {
                        Dx_vect[ns - 1] = Dx_vect[ns - 1] + shift_Dx_best;
                        Dy_vect[ns - 1] = Dy_vect[ns - 1] + shift_Dy_best;
                    }
                }
                else
                {
                    //Dx_vect[ns - 1] = old_Dx_vect[ns - 1];
                    //Dy_vect[ns - 1] = old_Dy_vect[ns - 1];
                }
            }

        }



    }
}