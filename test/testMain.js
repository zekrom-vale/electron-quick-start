// Tests main.js creation and clean up
const assert = require('assert/strict');
const test = require('node:test');
const {app, BrowserWindow} = require('electron')

console.log("Main Test")
console.log(`path: ${app.getAppPath()}`)
test("Main Test", async (contex) =>{
	test.after(app.exit)
	const main = require("../main.js")
	
	app.on('ready', main.createWindow)

	// Quit when all windows are closed.
	app.on('window-all-closed', function () {

	  console.log(now()+'::window-all-closed')
	  main.cleanUpApplication()
	  if(config.get("app.quitOnClose"))app.quit()
	})

	await new Promise((r,x)=>{setTimeout(r, 100)})
	
	// Hook into R process standard error
	var pass = true
    main.childProcess.stderr.on('data', data => {
		// Fix array of chars
		data = String.fromCharCode(...data)
		if(/halted/ig.test(data)){
			pass = false
			console.error("Failing: Halted")
			assert.fail(`R Halted: ${data}`)
		}
		if(/error/ig.test(data)){
			pass = false
			console.warn("Failing: Error")
			assert.fail(`R Had an error: ${data}`)
		}
		if(/warning/ig.test(data)){
			console.warn("Warning Detected")
		}
	})
	
	await new Promise((r,x)=>setTimeout(r, 7*1000))
	
	// Test reload
	main.mainWindow.reload()
	
	await new Promise((r,x)=>setTimeout(r, 7*1000))
	
	main.cleanUpApplication()
	// TODO: Add test to check if cleanup passed
	
	await new Promise((r,x)=>setTimeout(r, .5*1000))
	
	assert.ok(pass)
})