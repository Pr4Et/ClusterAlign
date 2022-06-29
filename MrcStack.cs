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
using System.IO;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;
using Emgu.CV;
using Emgu.Util;
using Emgu.CV.Structure;
using System.Drawing;
using System.Runtime;
using Emgu.CV.CvEnum;

namespace ClusterAlign
{
    public class MrcStack
    {
        public const int MRC_MODE_BYTE = 0;
        public const int MRC_MODE_SHORT = 1;
        public const int MRC_MODE_FLOAT = 2;
        public const int MRC_MODE_COMPLEX_SHORT = 3;
        public const int MRC_MODE_COMPLEX_FLOAT = 4;
        public const int MRC_MODE_USHORT = 6;
        public const int MRC_MODE_RGB = 16;
        public const int MRC_CFLOAT_REAL_IMAG = 20;         /* COMPLEX FLOAT mode */
public const int MRC_CFLOAT_AMP_RAD = 21;        /* COMPLEX FLOAT mode, but in amplitude and phase(  ) form */
        public const int MRC_CFLOAT_AMP_DEG = 22;        /* COMPLEX FLOAT mode, but in amplitude and phase( degree ) form */
        public const int MRC_LABEL_SIZE = 80;
        public const int MRC_NEXTRA = 16;
        public const int MRC_NLABELS = 10;
        public const int MRC_HEADER_SIZE = 1024;  /* Length of Header is 1024 Bytes. */
        public const int MRC_MAXCSIZE = 3;
        public const int headersize = 1024; //actually it is not accurate if nint!=0 or nreal!=0, see https://bio3d.colorado.edu/imod/doc/mrc_format.txt
        public static string headerType;
        public static bool imodstack;
        /* The header structure for MRC files, see: https://www.ccpem.ac.uk/mrc_format/mrc2014.php  */
        [StructLayout(LayoutKind.Explicit, Size = headersize, Pack = 1)]
        public unsafe struct MRCheader
        {
            [FieldOffset(0)]
            public int nx;         /*  # of Columns                  */
            [FieldOffset(4)]
            public int ny;         /*  # of Rows                     */
            [FieldOffset(8)]
            public int nz;         /*  # of Sections.                */
            [FieldOffset(12)]
            public int mode;       /*  given by #define MRC_MODE...  */
            [FieldOffset(16)]
            public int nxstart;    /*  Starting point of sub image.  */
            [FieldOffset(20)]
            public int nystart;
            [FieldOffset(24)]
            public int nzstart;
            [FieldOffset(28)]
            public int mx;         /* Grid size in x, y, and z       */
            [FieldOffset(32)]
            public int my;
            [FieldOffset(36)]
            public int mz;
            [FieldOffset(40)]
            public float xlen;       /* length of x element in um.     */
            [FieldOffset(44)]
            public float ylen;       /* get scale = xlen/nx ...        */
            [FieldOffset(48)]
            public float zlen;
            [FieldOffset(52)]
            public float alpha;      /* cell angles, ignore */
            [FieldOffset(56)]
            public float beta;
            [FieldOffset(60)]
            public float gamma;
            [FieldOffset(64)]
            public int mapc;       /* map coloumn 1=x,2=y,3=z.       */
            [FieldOffset(68)]
            public int mapr;       /* map row     1=x,2=y,3=z.       */
            [FieldOffset(72)]
            public int maps;       /* map section 1=x,2=y,3=z.       */
            [FieldOffset(76)]
            public float dmin;
            [FieldOffset(80)]
            public float dmax;
            [FieldOffset(84)]
            public float dmean;
            [FieldOffset(88)]
            public short ispg;       /* image type */
            [FieldOffset(90)]
            public short nsymbt;     /* space group number */
            /* 64 bytes */
            [FieldOffset(92)]
            public int next;
            [FieldOffset(96)]
            public short creatid;  /* Creator id, hvem = 1000, DeltaVision = -16224 */
            [FieldOffset(98)]
            public fixed byte blank[6];
            [FieldOffset(104)]
            public fixed byte extType[4]; //'SERI' for SerialEM, 'FEI1' for FEI, 'AGAR' for Agard
            [FieldOffset(108)]
            public int NVersion;
            [FieldOffset(112)]
            public fixed byte blank2[16];
            [FieldOffset(128)]
            public short nint;
            [FieldOffset(130)]
            public short nreal;
            [FieldOffset(132)]
            public short sub;
            [FieldOffset(134)]
            public short zfac;
            [FieldOffset(136)]
            public float min2;
            [FieldOffset(140)]
            public float max2;
            [FieldOffset(144)]
            public float min3;
            [FieldOffset(148)]
            public float max3;
            [FieldOffset(152)]
            public fixed byte imodStamp[4]; //'IMOD' stamp after creating fixed stack in IMOD
            [FieldOffset(156)]
            public float max4;
            [FieldOffset(160)]
            public short idtype;
            [FieldOffset(162)]
            public short lens;
            [FieldOffset(164)]
            public short nd1;     /* Devide by 100 to get float value. */
            [FieldOffset(166)]
            public short nd2;
            [FieldOffset(168)]
            public short vd1;
            [FieldOffset(170)]
            public short vd2;
            [FieldOffset(172)]
            public fixed float tiltangles[6];//=new float[6];  /* 0,1,2 = original:  3,4,5 = current */
            [FieldOffset(196)]
            public float xorg;
            [FieldOffset(200)]
            public float yorg;
            [FieldOffset(204)]
            public float zorg;
            [FieldOffset(208)]
            public fixed byte cmap[4];  
            [FieldOffset(212)]
            public fixed byte stamp[4]; //machine stamp: little-endian:0x44 0x44 0x00 0x00, big-endian: 0x11 0x11 0x00 0x00
            [FieldOffset(216)]
            public float rms;
            [FieldOffset(220)]
            public int nlabl;
            [FieldOffset(224)]
            public fixed byte labels1[80]; //in c++ char is the size of byte, was char  labels[10][80]
            [FieldOffset(304)]
            public fixed byte labels2[80];
            [FieldOffset(384)]
            public fixed byte labels3[80];
            [FieldOffset(464)]
            public fixed byte labels4[80];
            [FieldOffset(544)]
            public fixed byte labels5[80];
            [FieldOffset(624)]
            public fixed byte labels6[80];
            [FieldOffset(704)]
            public fixed byte labels7[80];
            [FieldOffset(784)]
            public fixed byte labels8[80];
            [FieldOffset(864)]
            public fixed byte labels9[80];
            [FieldOffset(944)]
            public fixed byte labels10[80];
        }
        [StructLayout(LayoutKind.Explicit,Size=headersize,Pack =1)]
        public unsafe struct headerAccess
        {
            [FieldOffset(0)]
            public MRCheader header;     // for field recognition
            [FieldOffset(0)]
            public fixed byte inbytes[headersize];
        }
 
        public unsafe headerAccess readheader;
        public FileStream infile,outfile;
        public string toshowfilename;
        
        public bool is_cashed;


 	        public int Size() {return readheader.header.nz;}

            public int Width() {return readheader.header.nx;}

            public int Height() {return readheader.header.ny;}

        public unsafe Mat[] FOpen(string filename) {
            infile = File.OpenRead(filename);
            toshowfilename = filename;
            byte[] arr = new byte[headersize];
            infile.Read(arr, 0, headersize); //read into header
            fixed (byte* byteptr= readheader.inbytes)
            { Marshal.Copy(arr, 0, new IntPtr(byteptr), headersize); } //Marshal.Copy(byte[] source, int startIndex, IntPtr destination, int length);
            PrintHeader();
            Mat[] slices=new Mat[readheader.header.nz];

            for (uint i = 0; i < readheader.header.nz; i++)
            {
                slices[i]= ReadImage(i);
            }
            infile.Close();
            return slices;
        }

        public unsafe void Fexport(string ref_filename, string out_filename, ref double[] Dx_vect, ref double[] Dy_vect, ref double report_phi)
        {
            if (!Program.isMRCfile)
            {
                Console.WriteLine("Aligned image file cannot be generated from dataset in TIF format.");
                return;
            }
            bool rotate2normal = (Settings4ClusterAlign2.Default.export_normalali || (Settings4ClusterAlign2.Default.isArbAngle && !out_filename.Contains("_jali."))) && (report_phi!=0);
            float phi = (float)report_phi;
            outfile = File.OpenWrite(out_filename);
            byte[] arr = new byte[headersize];
            infile = File.OpenRead(ref_filename);
            infile.Read(arr, 0, headersize); //read into header
            fixed (byte* byteptr = readheader.inbytes)
            { Marshal.Copy(arr, 0, new IntPtr(byteptr), headersize); } //Marshal.Copy(byte[] source, int startIndex, IntPtr destination, int length);
            int nx = readheader.header.nx;
            int ny = readheader.header.ny;
            infile.Close();
            infile = File.OpenRead(ref_filename);
            byte[] outarr = new byte[headersize + readheader.header.next];
            infile.Read(outarr, 0, headersize + readheader.header.next); //read into header
            outfile.Write(outarr, 0, headersize + readheader.header.next); //write header
            
            Mat slice;
            Mat warpMat = new Mat(2, 3, DepthType.Cv32F, 1);
            //CvInvoke.GetRotationMatrix2D(new PointF((float)0.5*ny,(float)0.5*nx), 0, 1, warpMat);//center, angle, scale
            PointF[] srcTri=new PointF[3];
            PointF[] dstTri=new PointF[3];

            for (uint i = 0; i < readheader.header.nz; i++)
            {
                slice = ReadImage(i);
                float x0 = (float)0.25 * nx + (float)Dx_vect[i];
                float y0 = (float)0.25 * ny + (float)Dy_vect[i];
                float x1 = (float)0.75 * nx + (float)Dx_vect[i];
                float y1 = (float)0.25 * ny + (float)Dy_vect[i];
                float x2 = (float)0.75 * nx + (float)Dx_vect[i];
                float y2 = (float)0.75 * ny + (float)Dy_vect[i];
                srcTri[0] = new PointF((float)0.25 * ny, (float)0.25 * nx);
                srcTri[1] = new PointF((float)0.25 * ny, (float)0.75 * nx);
                srcTri[2] = new PointF((float)0.75 * ny, (float)0.75 * nx);
                dstTri[0] = new PointF(y0, x0);
                dstTri[1] = new PointF(y1, x1);
                dstTri[2] = new PointF(y2, x2);
                if (rotate2normal)
                {
                    dstTri[0] = new PointF(0.5f * ny + MathF.Cos(-phi) * (y0 - 0.5f*ny) - MathF.Sin(-phi) * (x0 - 0.5f*nx), 0.5f*nx + MathF.Cos(-phi) * (x0 - 0.5f*nx) + MathF.Sin(-phi) * (y0 - 0.5f*ny));
                    dstTri[1] = new PointF(0.5f * ny + MathF.Cos(-phi) * (y1 - 0.5f*ny) - MathF.Sin(-phi) * (x1 - 0.5f*nx), 0.5f*nx + MathF.Cos(-phi) * (x1 - 0.5f*nx) + MathF.Sin(-phi) * (y1 - 0.5f*ny));
                    dstTri[2] = new PointF(0.5f * ny + MathF.Cos(-phi) * (y2 - 0.5f*ny) - MathF.Sin(-phi) * (x2 - 0.5f*nx), 0.5f*nx + MathF.Cos(-phi) * (x2 - 0.5f*nx) + MathF.Sin(-phi) * (y2 - 0.5f*ny));
                }
                Mat warp_mat = CvInvoke.GetAffineTransform(srcTri, dstTri);
                MCvScalar border_shade=CvInvoke.Mean(slice);
                CvInvoke.WarpAffine(slice, slice, warp_mat, slice.Size,Inter.Linear,Warp.Default, BorderType.Constant,border_shade);
                WriteImage(i, slice);
            }
            infile.Close();
            outfile.Close();
            //return true;
        }

 
        private void WriteImage(uint index, Mat slice)
        {
            int nx = readheader.header.nx;
            int ny = readheader.header.ny;

            switch (readheader.header.mode)
            {
                case MRC_MODE_BYTE:
                    {
                        int bufsize = readheader.header.nx * readheader.header.ny; 
                        byte[] tmpbuf = new byte[bufsize];
                        Array arr_slice = new byte[ny, nx]; //Nrows, Ncols = y,x
                        arr_slice = slice.GetData();  //from mat to array
                        byte temp;
                        for (int y = 0; y < ny; y++)
                        {
                            for (int x = 0; x < nx ; x++)
                            {
                                temp = Convert.ToByte(arr_slice.GetValue(x, y));
                                tmpbuf[x + y * nx] = temp;
                            }
                        }
                        outfile.Write(tmpbuf, 0, bufsize); //buffer, offset in buffer, size
                        break;
                    }
                case MRC_MODE_SHORT:
                    {
                        int bufsize = readheader.header.nx * readheader.header.ny * sizeof(short); //two bytes in int16
                        byte[] tmpbuf = new byte[bufsize];
                        Array arr_slice = new short[ny, nx]; //Nrows, Ncols = y,x
                        arr_slice = slice.GetData();  //from mat to array
                        byte[] temp = new byte[2];
                        for (int y = 0; y < ny; y++)
                        {
                            for (int x = 0; x < nx ; x++)
                            {
                                temp = BitConverter.GetBytes(Convert.ToInt16(arr_slice.GetValue(x, y)));
                                tmpbuf[x * sizeof(short) + y * nx * sizeof(short)] = temp[0];
                                tmpbuf[1 + x * sizeof(short) + y * nx * sizeof(short)] = temp[1];
                            }
                        }
                        outfile.Write(tmpbuf, 0, bufsize); //buffer, offset in buffer, size
                        break;
                    }
                case MRC_MODE_USHORT:
                    {
                        int bufsize = readheader.header.nx * readheader.header.ny * sizeof(short); //two bytes in uint16
                        byte[] tmpbuf = new byte[bufsize];
                        Array arr_slice = new ushort[ny, nx]; //Nrows, Ncols = y,x
                        arr_slice = slice.GetData();  //from mat to array
                        byte[] temp = new byte[2];
                        for (int y = 0; y < ny; y++)
                        {
                            for (int x = 0; x < nx; x++)
                            {
                                temp = BitConverter.GetBytes(Convert.ToUInt16(arr_slice.GetValue(x, y)));
                                tmpbuf[x * sizeof(short) + y * nx * sizeof(short)] = temp[0];
                                tmpbuf[1 + x * sizeof(short) + y * nx * sizeof(short)] = temp[1];
                            }
                        }
                        outfile.Write(tmpbuf, 0, bufsize); //buffer, offset in buffer, size
                        break;
                    }
                case MRC_MODE_FLOAT:
                    {
                        int bufsize = readheader.header.nx * readheader.header.ny * sizeof(float); //4 bytes in float
                        byte[] tmpbuf = new byte[bufsize];
                        Array arr_slice = new float[ny, nx]; //Nrows, Ncols = y,x
                        arr_slice = slice.GetData();  //from mat to array
                        byte[] temp = new byte[4];
                        for (int y = 0; y < ny; y++)
                        {
                            for (int x = 0; x < nx ; x++)
                            {
                                temp = BitConverter.GetBytes(Convert.ToSingle(arr_slice.GetValue(x, y)));
                                tmpbuf[x * sizeof(float) + y * nx * sizeof(float)] = temp[0];
                                tmpbuf[1 + x * sizeof(float) + y * nx * sizeof(float)] = temp[1];
                                tmpbuf[2 + x * sizeof(float) + y * nx * sizeof(float)] = temp[2];
                                tmpbuf[3 + x * sizeof(float) + y * nx * sizeof(float)] = temp[3];
                            }
                        }
                        outfile.Write(tmpbuf, 0, bufsize); //buffer, offset in buffer, size
                        break;
                    }
                default:
                    { break; }
            }


        }

        private Mat ReadImage(uint index)
        {
            int nx = readheader.header.nx;
            int ny = readheader.header.ny;
            Size size=new Size(readheader.header.nx, readheader.header.ny);

 
            IntPtr imptr= CvInvoke.cvCreateMat(nx, ny, DepthType.Cv32F);

            int bufsize = readheader.header.nx * readheader.header.ny;

            switch (readheader.header.mode)
            {
            case MRC_MODE_BYTE:
                {
                    
                    byte[] tmpbuf = new byte[bufsize];

                    for (int y = 0; y < ny; y++)
                    {
                    infile.Seek(headersize + readheader.header.next + (index * nx * ny + y * nx) * sizeof(byte), SeekOrigin.Begin); //AppDomainSetup the read pointer at location
                    infile.Read(tmpbuf, 0, nx * sizeof(byte)); //buffer, offset in buffer, size
                            for (int x = 0; x < nx; x++) {
                                CvInvoke.cvSetReal2D(imptr, x, y,tmpbuf[x] / 255.0f); //convert from byte to the opencv data already defined float
                            }
                    }
                break;
                 }
            case MRC_MODE_SHORT: //16 bits signed
                 {
                        byte[] tmpbuf = new byte[bufsize];
                        for (int y = 0; y < ny; y++)
                        {
                            infile.Seek(headersize + readheader.header.next + (index * nx * ny + y * nx) * sizeof(short), SeekOrigin.Begin); //AppDomainSetup the read pointer at location
                            infile.Read(tmpbuf, 0, nx * sizeof(short)); //buffer, offset in buffer, size
                            for (int x = 0; x < nx * sizeof(short); x += 2)
                            {
                                CvInvoke.cvSetReal2D(imptr, x/sizeof(short), y, (BitConverter.ToInt16(tmpbuf, x))); //two bytes to int16
                            }
                        }
                        break;
                 }
            case MRC_MODE_USHORT: //16 bits unsigned
                 {
                        byte[] tmpbuf = new byte[bufsize];
                        for (int y = 0; y < ny; y++)
                        {
                            infile.Seek(headersize + readheader.header.next + (index * nx * ny + y * nx) * sizeof(short), SeekOrigin.Begin); //AppDomainSetup the read pointer at location
                            infile.Read(tmpbuf, 0, nx * sizeof(short)); //buffer, offset in buffer, size
                            for (int x = 0; x < nx * sizeof(short); x += 2)
                            {
                                CvInvoke.cvSetReal2D(imptr, x/ sizeof(short), y, (BitConverter.ToUInt16(tmpbuf, x))); //two bytes to uint16
                            }
                        }
                        break;
                 }
            case MRC_MODE_FLOAT:
                 {
                        byte[] tmpbuf = new byte[bufsize];
                        for (int y = 0; y < ny; y++)
                        {
                            infile.Seek(headersize + readheader.header.next + (index * nx * ny + y * nx) * sizeof(float), SeekOrigin.Begin); //AppDomainSetup the read pointer at location
                            infile.Read(tmpbuf, 0, nx * sizeof(float)); //buffer, offset in buffer, size
                            for (int x = 0; x < nx * sizeof(float); x += sizeof(float))
                            {
                                CvInvoke.cvSetReal2D(imptr, x/ sizeof(float), y, (BitConverter.ToSingle(tmpbuf, x))); //4bytes bytes to float
                            }
                        }

                        break;
                 }
            default:
                 { break; }
             }

         return CvInvoke.CvArrToMat(imptr);
    }


        public unsafe void PrintHeader()
        {
            Console.WriteLine("");
            Console.WriteLine("mrc file information: "+ toshowfilename);
            Console.WriteLine("nx: " + readheader.header.nx.ToString());         /*  # of Columns                  */
            Console.WriteLine("ny: " + readheader.header.ny.ToString());         /*  # of Columns                  */
            Console.WriteLine("nz: " + readheader.header.nz.ToString());         /*  # of Columns                  */
            Console.WriteLine("mode: " + readheader.header.mode.ToString());         /*  MRC mode                 */
            //Console.WriteLine("nxstart: " + readheader.header.nxstart.ToString());         /*  Starting point of sub image.  */
            //Console.WriteLine("nystart: " + readheader.header.nystart.ToString());         /*                 */
            //Console.WriteLine("nzstart: " + readheader.header.nzstart.ToString());         /*                 */
            //Console.WriteLine("mx: " + readheader.header.mx.ToString());         /* Grid size in x, y, and z               */
            //Console.WriteLine("my: " + readheader.header.my.ToString());         /*   */
            //Console.WriteLine("mz: " + readheader.header.mz.ToString());         /*               */
            //Console.WriteLine("alpha: " + readheader.header.alpha.ToString());         /* cell angles          */
            //Console.WriteLine("beta: " + readheader.header.beta.ToString());         /*                */
            //Console.WriteLine("gamma: " + readheader.header.gamma.ToString());         /*                 */
            //Console.WriteLine("mapc: " + readheader.header.mapc.ToString());         /*  map coloumn 1=x,2=y,3=z.                  */
            //Console.WriteLine("mapr: " + readheader.header.mapr.ToString());         /*            */
            //Console.WriteLine("maps: " + readheader.header.maps.ToString());         /*                  */
            Console.WriteLine("dmin: " + readheader.header.dmin.ToString());         /*  # of Columns                  */
            Console.WriteLine("dmax: " + readheader.header.dmax.ToString());         /*  # of Columns                  */
            Console.WriteLine("dmean: " + readheader.header.dmean.ToString());         /*  # of Columns                  */
            //Console.WriteLine("ispg: " + readheader.header.ispg.ToString());         /*  #image type                 */
            //Console.WriteLine("nsymbt: " + readheader.header.nsymbt.ToString());         /*  space group number              */
            Console.WriteLine("MRC Version: " + readheader.header.NVersion.ToString());         /*  space group number              */
            string Result = "";
            string HexAlphabet = "0123456789ABCDEF";
            byte[] B=new byte[4];
            fixed (byte* byteptr = readheader.header.stamp)
            { Marshal.Copy(new IntPtr(byteptr), B, 0, 4); } //Marshal.Copy(IntPtr source, byte[] destination, int startIndex, int length);
            for (int n=0;n<4;n++)
            {
                Result+=HexAlphabet[(int)(B[n] >> 4)];
                Result += HexAlphabet[(int)(B[n] & 0xF)]+" ";
            }
            Console.WriteLine("machine-stamp: " + Result);         /*  space group number              */
            fixed (byte* byteptr = readheader.header.extType)
            { Marshal.Copy(new IntPtr(byteptr), B, 0, 4); } //Marshal.Copy(IntPtr source, byte[] destination, int startIndex, int length);
            headerType = Encoding.Default.GetString(B);
            Console.WriteLine("Header Type: " + headerType);
            if (readheader.header.nint>0 || readheader.header.nreal>0)
            {
                Console.WriteLine("Warning: Header extension is not supported. nint=" + readheader.header.nint.ToString()+"  nreal="+ readheader.header.nreal.ToString());
            }
            fixed (byte* byteptr = readheader.header.imodStamp)
            { Marshal.Copy(new IntPtr(byteptr), B, 0, 4); } //Marshal.Copy(IntPtr source, byte[] destination, int startIndex, int length);
            imodstack = (Encoding.Default.GetString(B)=="IMOD");
            if (!imodstack && headerType=="FEI1")
            {
                Console.WriteLine("### WARNING:  If you intend to use IMOD you first need to create fixed stack, since IMOD corrects the coordinate system of the FEI microscope. ###");
            }
        }

    }
}
