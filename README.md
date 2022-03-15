***This is the second revision to combine Electron with [R-Portable](https://sourceforge.net/projects/rportable/files/R-Portable/) (version 4.1.1) and [RStudio Shiny](https://cran.r-project.org/web/packages/shiny/index.html) package to deliver Electron Applications that run standalone with R and Shiny.***


# electron-quick-start

**Clone and run for a quick way to see Electron in action with R's Shiny.**

This is a minimal Electron and R application that expands on the [Quick Start Guide](https://electronjs.org/docs/tutorial/quick-start)
within the Electron documentation.  This allows any shiny application to be run portably without having R installed directly on the user's computer.

**Use this app along with the [Electron API Demos](https://electronjs.org/#get-started) app for API code examples to help you get started.**

## Files
### A basic Electron application needs just these files:
- `package.json` - Points to the app's main file and lists its details and dependencies.
- `main.js` - Starts the app and creates a browser window to render HTML. This is the app's **main process**.
- `index.html` - A web page to render. This is the app's **renderer process**.

You can learn more about each of these components within the [Quick Start Guide](https://electronjs.org/docs/tutorial/quick-start).

### R adds a few more files required, these are used and defined in `main.js`:
- `R-Portable-Win` - The windows R files (Just take them from ProgramFiles)
- `R-Portable-Mac` - The mac R files (Just take them from Applications and run `R-Portable-Mac/bin/R` to untangle it)
- `R-Portable-*/library` - The location to protably install R packages
- `app.R` - The shiny application to start
- `cc.ico` - The icon used for the application made
- Other files to be coppied to `resources\app`


## To Use
### Installing
To clone and run this repository you'll need [Git](https://git-scm.com) and [Node.js](https://nodejs.org/en/download/) (which comes with [npm](http://npmjs.com)) installed on your computer. From your command line:

```bash
# Clone this repository
git clone https://github.com/zekrom-vale/electron-quick-start
# Install Electron Packager (if first time)
npm install electron-packager -g 
# Go into the repository
cd electron-quick-start
# Install dependencies
npm install
```

Then install any pakages with R (Not with the portable version as that does not resolve dependencies).  Copy those installed files to the relevent library `R-Portable-*/library`.  Do not remove the existing library files as it breaks the base packages.

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

#### Build the Executable/App for macOS
```bash
cd electron-quick-start
npm run package-mac
```

## Notes
 - To see the console that Electron prints out you can simply run the app from the command line.
 - If you're using Linux Bash for Windows, [see this guide](https://www.howtogeek.com/261575/how-to-run-graphical-linux-desktop-applications-from-windows-10s-bash-shell/) or use `node` from the command prompt.
 - macOS and linux support is experamental and may need some modification

## Forked from
 - [COVAIL/electron-quick-start](https://github.com/COVAIL/electron-quick-start) - The older version that is not beeing updated

## Resources for Learning Electron


- [electron/electron-quick-start](https://github.com/electron/electron-quick-start) - a very basic starter Electron app
- [electron/simple-samples](https://github.com/electron/simple-samples) - small applications with ideas for taking them further
- [electron/electron-api-demos](https://github.com/electron/electron-api-demos) - an Electron app that teaches you how to use Electron
- [hokein/electron-sample-apps](https://github.com/hokein/electron-sample-apps) - small demo apps for the various Electron APIs

## License

[CC0 1.0 (Public Domain)](LICENSE.md)
