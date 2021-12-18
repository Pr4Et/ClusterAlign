// This code for establishing Avalonia GUI is adapted from Wiesław Šoltés under MIT license  
using Avalonia;
using Avalonia.Markup.Xaml;
using Avalonia.Controls.ApplicationLifetimes;
using Avalonia.Logging.Serilog;


namespace ClusterAlign
{
    /// <summary>
    /// Main application.
    /// </summary>
    public class App : Application
    {
        static IClassicDesktopStyleApplicationLifetime mydesktopLifetime;
        /// <summary>
        /// Program entry point.
        /// </summary>
        /// <param name="args">The program arguments.</param>
        static void Main(string[] args)
        {
            BuildAvaloniaApp().StartWithClassicDesktopLifetime(args);
        }

        /// <summary>
        /// Build Avalonia app.
        /// </summary>
        /// <returns>The Avalonia app builder.</returns>
        public static AppBuilder BuildAvaloniaApp()
            => AppBuilder.Configure<App>()
                         .UsePlatformDetect()
                         .LogToDebug();

        /// <inheritdoc/>
        public override void Initialize()
        {
            AvaloniaXamlLoader.Load(this);
        }

        /// <inheritdoc/>
        public override void OnFrameworkInitializationCompleted()
        {
            if (ApplicationLifetime is IClassicDesktopStyleApplicationLifetime desktopLifetime)
            {
                desktopLifetime.MainWindow = new Window1();
                mydesktopLifetime = desktopLifetime;
            }
            else if (ApplicationLifetime is ISingleViewApplicationLifetime singleViewLifetime)
            {
                //singleViewLifetime.MainView = new Window1();
            }
            base.OnFrameworkInitializationCompleted();
        }
        public static void Hidewin()
        {
            mydesktopLifetime.MainWindow.Hide();
 
        }
    }
}