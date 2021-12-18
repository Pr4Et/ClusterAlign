//Shahar Seifer 2020
using Avalonia;
using Avalonia.Controls;
using Avalonia.Markup.Xaml;
using Avalonia.Interactivity;
using System.Threading.Tasks;
using System.IO;

namespace ClusterAlign
{
    public class Window1 : Window
    {
        public static string loadAttentionFile = "";
        public static string loadfidFile = "";
        public static int AttentionNumMarkedfids;

        private TextBox TtextInputPath;
        private Button btnBrowsePath;
        private Button btnloadAttention;
        private Button btnloadfid;
        private TextBox txtuseAttention;
        private TextBox txtMarkedfids;
        private NumericUpDown NumMarkedfids;
        private Button btnRun;
        private TextBox TDatafile;
        private TextBox TTiltFile;
        private ComboBox CombxisRotation;
        private NumericUpDown NumClusterSize;
        private NumericUpDown NumNcluster;
        private NumericUpDown Numfidsize;
        private NumericUpDown NumNfidmax;
        private NumericUpDown Numncenter;
        private NumericUpDown NumNminimumTracked;
        private ComboBox CombForceFill;
        private ComboBox Combcoswindow;
        private NumericUpDown NumPreAlignmentTol;
        private ComboBox Combfiducialbright;
        private NumericUpDown NumTolFidCenter;
        private NumericUpDown NumTolFidSize;
        public Window1()
        {
            this.InitializeComponent();
            #if DEBUG
                    this.AttachDevTools();
            #endif

            TtextInputPath = this.FindControl<TextBox>("textInputPath");
            btnBrowsePath = this.FindControl<Button>("btnBrowsePath");
            btnloadAttention = this.FindControl<Button>("loadAttention");
            btnloadfid = this.FindControl<Button>("loadfid");
            txtuseAttention = this.FindControl<TextBox>("useAttention");
            NumMarkedfids = this.FindControl<NumericUpDown>("NumMarkedfids");
            txtMarkedfids = this.FindControl<TextBox>("Markedfids"); ;
            TDatafile = this.FindControl<TextBox>("Datafile"); 
            TTiltFile = this.FindControl<TextBox>("TiltFile");
            CombxisRotation = this.FindControl<ComboBox>("xisRotation");
            NumClusterSize = this.FindControl<NumericUpDown>("ClusterSize");
            NumNcluster = this.FindControl<NumericUpDown>("Ncluster");
            Numfidsize = this.FindControl<NumericUpDown>("fidsize");
            NumNfidmax = this.FindControl<NumericUpDown>("NFidMax");
            Numncenter = this.FindControl<NumericUpDown>("ncenter");
            NumNminimumTracked = this.FindControl<NumericUpDown>("NminimumTracked");
            CombForceFill = this.FindControl<ComboBox>("ForceFill");
            Combcoswindow = this.FindControl<ComboBox>("coswindow");
            NumPreAlignmentTol = this.FindControl<NumericUpDown>("PreAlignmentTol");
            Combfiducialbright = this.FindControl<ComboBox>("fiducialbright");
            NumTolFidCenter = this.FindControl<NumericUpDown>("TolFidCenter");
            NumTolFidSize = this.FindControl<NumericUpDown>("TolFidSize");

            TtextInputPath.Text = ClusterAlign.Settings4ClusterAlign2.Default.Path;
            TDatafile.Text= ClusterAlign.Settings4ClusterAlign2.Default.DataFileName;
            TTiltFile.Text= ClusterAlign.Settings4ClusterAlign2.Default.TiltFileName;
            CombxisRotation.SelectedIndex = (ClusterAlign.Settings4ClusterAlign2.Default.xisRotation ? 1:0);
            NumClusterSize.Value = ClusterAlign.Settings4ClusterAlign2.Default.cluster_size;
            NumNcluster.Value = ClusterAlign.Settings4ClusterAlign2.Default.Ncluster;
            Numfidsize.Value = ClusterAlign.Settings4ClusterAlign2.Default.fidsize;
            NumNfidmax.Value = ClusterAlign.Settings4ClusterAlign2.Default.NfidMax;
            Numncenter.Value = ClusterAlign.Settings4ClusterAlign2.Default.ncenter;
            NumNminimumTracked.Value = ClusterAlign.Settings4ClusterAlign2.Default.N_minimum_tracked_fiducials;
            CombForceFill.SelectedIndex = (ClusterAlign.Settings4ClusterAlign2.Default.ForceFill ? 1 : 0);
            Combcoswindow.SelectedIndex = (ClusterAlign.Settings4ClusterAlign2.Default.coswindow ? 1 : 0);
            NumPreAlignmentTol.Value = ClusterAlign.Settings4ClusterAlign2.Default.PreAlignmentTol;
            NumTolFidSize.Value = ClusterAlign.Settings4ClusterAlign2.Default.TolFidSize;
            NumTolFidCenter.Value = ClusterAlign.Settings4ClusterAlign2.Default.TolFidCenter;
            Combfiducialbright.SelectedIndex = (ClusterAlign.Settings4ClusterAlign2.Default.fiducials_bright ? 1 : 0);
            txtuseAttention.Text ="Optional tif stack: N/A";
            txtMarkedfids.IsVisible = false;
            NumMarkedfids.Value = 0;
            NumMarkedfids.IsVisible = false;
            btnRun = this.FindControl<Button>("btnRun");
            btnBrowsePath.Click += async (sender, e) => await GetPath();
            btnloadAttention.Click += async (sender, e) => await GetAttentionFile();
            btnloadfid.Click += async (sender, e) => await GetfidFile();
            btnRun.Click += BtnRun_Click;

        }

        private void BtnRun_Click(object sender, RoutedEventArgs e)
        {
            AttentionNumMarkedfids = (int)NumMarkedfids.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.Path= TtextInputPath.Text;
            ClusterAlign.Settings4ClusterAlign2.Default.DataFileName= TDatafile.Text;
            ClusterAlign.Settings4ClusterAlign2.Default.TiltFileName= TTiltFile.Text;
            ClusterAlign.Settings4ClusterAlign2.Default.xisRotation = (CombxisRotation.SelectedIndex==1 ? true : false);
            ClusterAlign.Settings4ClusterAlign2.Default.cluster_size= (int)NumClusterSize.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.Ncluster= (int)NumNcluster.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.NfidMax =(int)NumNfidmax.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.fidsize= (float)Numfidsize.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.ncenter= (int)Numncenter.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.N_minimum_tracked_fiducials= (int)NumNminimumTracked.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.ForceFill = (CombForceFill.SelectedIndex==1 ? true : false);
            ClusterAlign.Settings4ClusterAlign2.Default.coswindow = (Combcoswindow.SelectedIndex==1 ? true : false);
            ClusterAlign.Settings4ClusterAlign2.Default.PreAlignmentTol = NumPreAlignmentTol.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.fiducials_bright = (Combfiducialbright.SelectedIndex==1 ? true : false);
            ClusterAlign.Settings4ClusterAlign2.Default.TolFidSize = (int)NumTolFidSize.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.TolFidCenter = (int)NumTolFidCenter.Value;
            ClusterAlign.Settings4ClusterAlign2.Default.Save();
            ClusterAlign.App.Hidewin();
            Program.MyMain();
            

        }

        private void InitializeComponent()
        {
            AvaloniaXamlLoader.Load(this);
        }
        //private void  DataFileButtonClick(object sender, RoutedEventArgs e)
        //{
        //    GetPath();
        //}
        private async Task GetPath()
        {
            //var dlg = new OpenFolderDialog();
            var dlg = new OpenFileDialog();
            string text = TtextInputPath.Text;
            if (text.Length > 0)
            {
                dlg.Directory = TtextInputPath.Text;
            }

            var result = await dlg.ShowAsync(this);
            if (!string.IsNullOrWhiteSpace(result[0]))
            {
                TtextInputPath.Text = Path.GetDirectoryName(@result[0]);
                TDatafile.Text=Path.GetFileName(@result[0]);
                TTiltFile.Text= Path.GetFileNameWithoutExtension(@result[0])+".rawtlt";
            }
        }

        private async Task GetAttentionFile()
        {
            var dlg = new OpenFileDialog();
            dlg.Title="Choose tif stack file in which important fiducials are painted over with color 0 paintbrush";
            var result = await dlg.ShowAsync(this);
            if (!string.IsNullOrWhiteSpace(result[0]))
            {
                loadAttentionFile = result[0];
                txtuseAttention.Text = "Optional tif stack: loaded";
                txtMarkedfids.IsVisible = true;
                //NumMarkedfids.IsVisible = true;

            }
        }

        private async Task GetfidFile()
        {
            var dlg = new OpenFileDialog();
            dlg.Title = "Choose *.fid.txt file with partial table of fiducials";
            var result = await dlg.ShowAsync(this);
            if (!string.IsNullOrWhiteSpace(result[0]))
            {
                loadfidFile = result[0];
            }
        }

        //private void DataFilePointerEnter(object sender, PointerEventArgs e)
        //{
        //    //do something when pointer enters
        // }

    }
}
