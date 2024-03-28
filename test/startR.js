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
			`(1221*124-786/2)**3;log(25890123570891234)`
		]
	)
	var i = 0;
	childProcess.stdout.on('data', data =>{
		switch(i++){
			case 3:
				assert.match(`${data}`, /3\.443703e\+15/i)
				break
			case 4:
				assert.match(`${data}`, /37\.79264/i)
		}
	})
	childProcess.stderr.on('data', data => {
		assert.fail(data)
	})
	await new Promise((r,x)=>{
		assert.ok(true)
		childProcess.on("close", r)
	})
	console.log("end")
})
