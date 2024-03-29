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
	var p = new Promise((r,x)=>{
		setTimeout(r, 7*1000)
	})
	main.childProcess.stdout.on('data', data => console.log(`Rout: ${data}`))
	var pass = true
    main.childProcess.stderr.on('data', data => {
		// Fix array of chars
		data = String.fromCharCode(...data)
		if(/halted/ig.test(data)){
			pass = false
			console.warn("Failing")
			assert.fail(`R Halted: ${data}`)
		}
		if(/error/ig.test(data)){
			pass = false
			console.warn("Failing")
			assert.fail(`R Had an error: ${data}`)
		}
	})
	await p
	
	main.cleanUpApplication()
	// TODO: Add test to check if cleanup passed
	assert.ok(pass)
})