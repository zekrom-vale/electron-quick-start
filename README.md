***This is the third revision to combine Electron with [R-Portable](https://sourceforge.net/projects/rportable/files/R-Portable/) (version 4.1.1) and [RStudio Shiny](https://cran.r-project.org/web/packages/shiny/index.html) package to deliver Electron Applications that run standalone with R and Shiny.***


# shinyElectron [![.github/workflows/buildTest.yml](https://github.com/zekrom-vale/shinyElectron/actions/workflows/buildTest.yml/badge.svg)](https://github.com/zekrom-vale/shinyElectron/actions/workflows/buildTest.yml)

**Clone and run for a quick way to see Electron in action with R's Shiny.**

This is a minimal Electron and R application that expands on the [Quick Start Guide](https://electronjs.org/docs/tutorial/quick-start)
within the Electron documentation.  This allows any shiny application to be run portably without having R installed directly on the user's computer.

**Use this app along with the [Electron API Demos](https://electronjs.org/#get-started) app for API code examples to help you get started.**

## Files
### A basic Electron application needs just these files:
- `package.json` - Points to the app's main file and lists its details and dependencies.
- `main.js` - Starts the app and creates a browser window to render HTML. This is the app's **main process**.

You can learn more about each of these components within the [Quick Start Guide](https://electronjs.org/docs/tutorial/quick-start).

### R adds a few more files required, these are used and defined in `main.js`:
- `R-Portable-Win` - The windows R files (Just take them from ProgramFiles)
- `R-Portable-Mac` - The mac R files (Just take them from Applications and run `R-Portable-Mac/bin/R` to untangle it)
- `R-Portable-Linux` - The linux R files (Build folowing [
Painless-R-compilation-...-Ubuntu](https://github.com/Jiefei-Wang/Painless-R-compilation-and-installation-on-Ubuntu) up to `make check` and compile packages now with `install.packages()` or `devtools`)
- `R-Portable-*/library` - The location to protably install R packages find these files with `.libPaths()` to find the library paths
- `app.R` - The shiny application to start
- `cc.ico` - The icon used for the application made
- `config` - Contains configuration files for building and runing the app
- `scripts` - Contains the instalation scripts used to build your application
- `DESCRIPTION` - Required for testing action on GitHub, can be safely removed
- Other files to be coppied to `resources\app`


## To Use
### Installing
To clone and run this repository you'll need [Git](https://git-scm.com) and [Node.js](https://nodejs.org/en/download/) (which comes with [npm](http://npmjs.com)) installed on your computer. From your command line:

```bash
# Clone this repository
git clone https://github.com/zekrom-vale/shinyElectron/tree/stable
# Install Electron Packager (if first time)
npm install -g --save-dev @electron/packager
# Go into the repository
cd electron-quick-start
# Install dependencies
npm install
```

Then install any pakages with R (Not with the portable version as that does not resolve dependencies).  Copy those installed files to the relevent library `R-Portable-*/library`.  Do not remove the existing library files as it breaks the base packages.

### Configure
The configuration files are included in the `config` folder, if you want to use aother file format feel free as long as [node-config](https://github.com/node-config/node-config/wiki/Configuration-Files#file-formats) supports it.

`default.yaml` The default options to use, see [node-config](https://github.com/node-config/node-config/wiki/Configuration-Files#file-load-order) to learn more on loading order.
```yaml
R:
    url: "http://127.0.0.1:" # The url of where shiny is hosted most likely the loop back at 127.0.0.1
    port: 9191 # The url of where shiny is hosted `${url}${port}`
    kill: false # Should R be killed on exit?
    app: app.R # The R script to run with shiny
    path:
        fixHome: true # Should the app fix R's `R_HOME_DIR` in `R-Portable-*/R` ignored if isPortable is false
        isPortable: true # should R be run as portable?
window:
    delay: 2000 # How long to wait to show the window
    poll: 1000 # How long to wait to try to conect to the URL again after failing
    loading:     # A URL to load for the loading bar
        path: loading.html # The URL or file to load
        isURL: false # Is path a file or a URL?
        config: # Loading window settings
            show: false # Show the window?  Using true may break things
            width: 1200 # Width of the window
            height: 1000 # Height of the window
            title: My app # Title of the loading window
    config: # Main window settings
        show: false # Show the window?  Using true may break things
        width: 1200 # Width of the window
        height: 1000 # Height of the window
        title: My app # Title of the window
        webPreferences: 
            nodeIntegration: false # Should node be interated into JS?  Not implimented yet
    dev: false # Load developer tools?
    fullReload: true # Reload the entire R sesion?
    # Recomended this to be true and R.kill be false
    # Use the folowing in R shiny server:
    # onSessionEnded(function(){
    #     quit(save = "no")
    #  })
app:
    CompanyName: None
    FileDescription: CE
    ProductName: Shiny Electron App
    out: ElectronShinyApp # Where to build the files to see below for more
    name: electron-quick-start # used in the specific arch and platform build under `${out}/${name}-${platform}-${arch}`
    icon: assets/icons/png/1024x1024.png # Icon of the run file
    quitOnClose: true # When all windows are closed is the application termniated?
```

`linux.yaml` Linux specific options, anything here will overwite `default.yaml`
```yaml
R:
    kill: pkill -9 "R" # The command to kill R
    path: ./R-Portable-Linux/bin/R # The path to R
    home: ./R-Portable-Linux # The home dir of R
```

`darwin.yaml` MacOS specific options
```yaml
R:
    kill: false # pkill -9 "R" # The command to kill R
    path: # The path to R
        portable: R-Portable-Mac/bin/R
        home: R-Portable-Mac # The home dir of R
        local: /Library/Frameworks/lib/R/bin/R 
app:
    quitOnClose: false # On macOS it's common to re-create a window in the app when the
                       # dock icon is clicked and there are no other windows open.
```

### Compile/Run
Then you can compile or run the app

#### Run the app
```bash
npm start
```

#### Build the Executable/App for Windows
```bash
cd electron-quick-start
npm run package-win
```

#### Build the Executable/App for Linux
```bash
cd electron-quick-start
npm run package-linux
```

#### Build the Executable/App for macOS
```bash
cd electron-quick-start
npm run package-mac
```

## Notes
 - To see the console that Electron prints out you can simply run the app from the command line.
 - If you're using Linux Bash for Windows, [see this guide](https://www.howtogeek.com/261575/how-to-run-graphical-linux-desktop-applications-from-windows-10s-bash-shell/) or use `node` from the command prompt.
 - macOS support is experamental and may need some modification

## Forked from
 - [COVAIL/electron-quick-start](https://github.com/COVAIL/electron-quick-start) - The older version that is not beeing updated

## Resources for Learning Electron


- [electron/electron-quick-start](https://github.com/electron/electron-quick-start) - a very basic starter Electron app
- [electron/simple-samples](https://github.com/electron/simple-samples) - small applications with ideas for taking them further
- [electron/electron-api-demos](https://github.com/electron/electron-api-demos) - an Electron app that teaches you how to use Electron
- [hokein/electron-sample-apps](https://github.com/hokein/electron-sample-apps) - small demo apps for the various Electron APIs

## License

[CC0 1.0 (Public Domain)](LICENSE.md)
