// Modules to control application life and create native browser window
const path = require('path')
const child = require('child_process')
process.env.NODE_ENV = process.platform
const config = require("config")
const assert = require('assert/strict');
const test = require('node:test');

const execPath = path.normalize(config.get("R.path.local"))

console.log("R Test")
test( "R Tests", async (contex) =>{
	var childProcess = child.spawn(
		execPath, [
			"-e",
			`1+1`
		]
	)
	childProcess.stdout.on('data', data => console.log(`Rout: ${data}`))
	childProcess.stderr.on('data', data => console.warn(`Rerr: ${data}`))
	await new Promise((r,x)=>{
		childProcess.on("close", r)
	})
	console.log("end")
})
