// Modules to control application life and create native browser window
const {app, BrowserWindow} = require('electron')
const path = require('path')
const shell = require('child_process').execSync
const to = require('await-to-js').default
const url = require('url')
const child = require('child_process')
process.env.NODE_ENV = process.platform
const config = require("config")

// Error
function error(msg, name=""){
  console.warn(`${name}: ${msg}`);
}

var appPath = path.normalize(config.get("R.app"))
	// This may not be required
	if(!path.isAbsolute(appPath))appPath=path.join(app.getAppPath(), appPath)
var port = parseInt(config.get("R.port"))
	if(is.NaN(port)){
		error(`R.port is not a Number change config settings got ${config.get("R.port")}, using port 9191`, "Port Warning")
		port = 9191
	}
const execPortable = !!config.get("R.path.isPortable")
const execPath = path.normalize(
		execPortable?
			path.join(app.getAppPath(), config.get("R.path.portable")):
			config.get("R.path.local")
	)

console.log(process.env)


//Platform switch
const MACOS = "darwin"
const WINDOWS = "win32"
const LINUX = "linux"
{
	// let and const here will be discarded after this block, var will be kept
	let _i=''
	switch(process.platform){
		case WINDOWS:
			break
		case MACOS:
			// Adapt to support MacOS sed
			_i=' ""'
			//Fall through as mac is linux like
		case LINUX:
			// Fix issue with R Home path by overriding the R sh script with the correct value
			if(execPortable && config.get("R.path.fixHome")){
				let home=path.join(app.getAppPath(), config.get("R.path.home"))
				// Must use ! as / is an issue with paths
				shell(`sed -i${_i} 's!R_HOME_DIR=.*$!R_HOME_DIR="${home}"!' ${execPath}`)
			}
			break
		default:
			error(`Not on windows, linux, or macos. Got ${process.platform}`, "Platform Error")
	}
}

// Looks like this is resolved
// Due to an issue with shiny, the port needs to be set via options and not passed to the runApp function
// https://github.com/rstudio/shiny/issues/1942
var childProcess = null
function sartR(){
    // Need to steralize R.port and appPath
    // Try to cast R.port to int
    childProcess = child.spawn(
    	execPath, [
    		"-e",
    		`options(shiny.port=${port});
shiny::runApp(file.path('${appPath}'))`
    	]
    )
    childProcess.stdout.on('data', data => console.log(`stdout:${data}`))
    childProcess.stderr.on('data', data => console.warn(`stderr:${data}`))
}

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
// To be updated later in createWindow
var mainWindow = null
async function createWindow(){
	sartR()
	
	console.log(now()+'create-window')
	let loading = new BrowserWindow(config.get("window.loading.config"))
    // May need to staralize the URL
    if(config.get("window.loading.isURL"))
    	loading.loadURL(config.get("window.loading.path"))
    else
    	loading.loadFile(path.normalize(config.get("window.loading.path")))
    
    if(config.get("window.dev"))loading.toggleDevTools()
	
	// loading.once('show', ...) only runs once
	
	////////////////////////////////////////////////////////////////////////////////
	// Show Loading
	////////////////////////////////////////////////////////////////////////////////
	
	{
		let p = new Promise((r, x)=>loading.once('show', r))
		loading.show()
		await p
	}
	
	mainWindow = new BrowserWindow(config.get("window.config"))

	////////////////////////////////////////////////////////////////////////////////
	// Connect to Shiny
	////////////////////////////////////////////////////////////////////////////////
	
	{
		let p = new Promise(
				(r, x)=>mainWindow.webContents.on('dom-ready',
					()=>setTimeout(r, config.get("window.delay"))
				)
			)
		
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
		// Open the DevTools if set in config
		if(config.get("window.dev"))mainWindow.webContents.openDevTools()
		
		////////////////////////////////////////////////////////////////////////////////
		// Shiny Ready
		////////////////////////////////////////////////////////////////////////////////
		
		await p
	}
	
	console.log(now()+'::mainWindow loaded')
	mainWindow.show()
	//Quit the loading page
	loading.hide()
	loading.close()
	
	////////////////////////////////////////////////////////////////////////////////
	// Window Closed or Reloaded
	////////////////////////////////////////////////////////////////////////////////
	// This may not work well
	if(config.get("window.fullReload") && config.get("R.kill"))mainWindow.webContents.on('beforeunload', function(){
		mainWindow.loadURL(config.get("window.loading"))
		setTimeout(()=>{
			console.log(now()+'::window-reload')
			console.log("==================================================================")
			console.log("Disposing prior window and R session, starting loading please wait")
			cleanUpApplication(false)
			createWindow()
		}, 1000)
	})
	// onSessionEnded(function(){
    // 		quit(save = "no")
  	// })
	if(config.get("window.fullReload") && !config.get("R.kill"))childProcess.on('close', function(){
		console.log(now()+'::R-close')
		console.log("==================================================================")
		console.log("Disposing prior window and R session, starting loading please wait")
		cleanUpApplication(false)
		createWindow()
	})
}

function now(){
	return new Date().toISOString()
}

function cleanUpApplication(quit=true){
	if(childProcess && config.get("R.kill")){
		// We will see if R gracefully exits
		childProcess.exit()
		// childProcess.kill()
		// We are not able to kill the program, not allowed 
		// child.execSync(kill)
	}
	if(quit){
		app.quit()
		// This will only exit when all async processes are finished
		// I could do process.kill(pid)
		process.exit()
	}
	else{
		mainWindow.hide()
		mainWindow.close()
	}
}
// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', createWindow)
app.on('activate', function () {
  // On macOS it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if(mainWindow === null) createWindow()
})

// Quit when all windows are closed.
app.on('window-all-closed', function () {

  console.log(now()+'::window-all-closed')
  cleanUpApplication()
	if(config.get("app.quitOnClose"))app.quit()
})


// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.

