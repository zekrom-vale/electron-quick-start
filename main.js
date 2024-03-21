// Modules to control application life and create native browser window
const {app, BrowserWindow} = require('electron')
const path = require('path')
const shell = require('child_process').execSync
const to = require('await-to-js').default
const url = require('url')
const child = require('child_process')

process.env.NODE_ENV = process.platform
const config = require("config")


const port = config.get("R.port")
const MACOS = "darwin"
const WINDOWS = "win32"
const LINUX = "linux"

if(process.platform == WINDOWS){
  appPath = appPath.replace(/\\/g, "\\\\");
}
else if(process.platform != LINUX && process.platform != MACOS) {
  console.log("not on windows or macos?")
  throw new Error("not on windows or macos?")
}

var appPath = config.get("R.app")
if(!path.isAbsolute(appPath)){
	appPath=path.join(app.getAppPath(), appPath)
}

var execPath
var execAbs
if(config.get("R.path.isPortable")){
	execAbs = true
	execPath = config.get("R.path.portable")
}
else{
	execAbs = false
	execPath = path.join(app.getAppPath(), config.get("R.path.local"))
}

console.log(process.env)

// Fix issue with R Home path
if(!execAbs && config.get("R.path.fixHome")){
	if(process.platform == LINUX){
		let home=path.join(app.getAppPath(), config.get("R.path.home") )
		shell(`sed -i 's!R_HOME_DIR=.*$!R_HOME_DIR="${home}"!' ${execPath}`)
	}
	else if(process.platform == MACOS){
		let home=path.join(app.getAppPath(), config.get("R.path.home") )
		shell(`sed -i "" 's!R_HOME_DIR=.*$!R_HOME_DIR="${home}"!' ${execPath}`)
	}
}

// Due to an issue with shiny, the port needs to be set via options and not passed to the runApp function
// https://github.com/rstudio/shiny/issues/1942
const childProcess = child.spawn(execPath, ["-e", `options(shiny.port=${port}); shiny::runApp(file.path('${appPath}'))`])
childProcess.stdout.on('data', (data) => {
  console.log(`stdout:${data}`)
})
childProcess.stderr.on('data', (data) => {
  console.log(`stderr:${data}`)
})

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow

function createWindow () {
  // Create the browser window.
  //  mainWindow = new BrowserWindow({webPreferences:{nodeIntegration:false},width: 800, height: 600})
  //  console.log(process.cwd())
  console.log('create-window')


    let loading = new BrowserWindow({show: false, frame: false})
    //let loading = new BrowserWindow()
    console.log(new Date().toISOString()+'::showing loading');
    loading.loadURL(config.get("window.loading"))
    ///loading.toggleDevTools()

    loading.once('show', async function(){
      console.log(new Date().toISOString()+'::show loading')
      mainWindow = new BrowserWindow(config.get("window.config"))
      mainWindow.webContents.once('dom-ready', () => {
        console.log(new Date().toISOString()+'::mainWindow loaded')
        setTimeout( () => {
          mainWindow.show()
          if(process.platform==MACOS){
            mainWindow.reload()
          }
          loading.hide()
          loading.close()
        }, config.get("window.delay"))
      })
      
	// Try loading the URL
	// Since Shiny is not instantly loaded, it needs to try to connect until it loads
	// Whithout this it likely will load a blank page until the user manualy reloads it
	{
	  let poll = config.get("window.poll")
	  while(true){
		  let [err, r] = await to(mainWindow.loadURL(config.get("R.url")+port))
		  //On failure it creates an async delay it waits for and continues the loop
		  if(err) await new Promise(r => setTimeout(r, poll))
		  else break
	  }
	}
      mainWindow.webContents.on('did-finish-load', function() {
        console.log(new Date().toISOString()+'::did-finish-load')
      });

      mainWindow.webContents.on('did-start-load', function() {
        console.log(new Date().toISOString()+'::did-start-load')
      });

      mainWindow.webContents.on('did-stop-load', function() {
        console.log(new Date().toISOString()+'::did-stop-load')
      });
      mainWindow.webContents.on('dom-ready', function() {
        console.log(new Date().toISOString()+'::dom-ready')
      });

      // Open the DevTools if set in config
      if(config.get("window.dev"))mainWindow.webContents.openDevTools()

      // Emitted when the window is closed.
      mainWindow.on('closed', function () {
        console.log(new Date().toISOString()+'::mainWindow.closed()')
        cleanUpApplication()
      })
    })

    loading.show()

}


function cleanUpApplication(){

  app.quit()
  let kill = config.get("R.kill")
  if(childProcess && kill){
    childProcess.kill()
    child.execSync(kill)      
  }
}
// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', createWindow)

// Quit when all windows are closed.
app.on('window-all-closed', function () {

  console.log('EVENT::window-all-closed')
  // On OS X it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  cleanUpApplication()
  // On macOS it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  //if (process.platform !== 'darwin') app.quit()

})

app.on('activate', function () {
  // On macOS it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (mainWindow === null) createWindow()
})

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.
