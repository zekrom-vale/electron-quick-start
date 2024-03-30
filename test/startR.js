// Tests the R part of the code
const path = require('path')
const child = require('child_process')
process.env.NODE_ENV = process.platform
const config = require("config")
const assert = require('assert/strict');
const {it, describe, after} = require('node:test');

const execPath = path.normalize(config.get("R.path.local"))

describe( "R Tests", async () =>{
	after(process.exit)
	
	var childProcess = child.spawn(
		execPath, [
			"-e",
			`(1221*124-786/2)**3;log(25890123570891234)`
		]
	)
	var i = 0;
	// Case may be shifted
	// Does not get to case 5: it(...)
	childProcess.stdout.on('data', data =>{
		switch(i++){
			case 4:
				it("Math 1", ()=>assert.match(`${data}`, /3\.443703\d*e\+15/i))
				break
			case 5:
				it("Math 2", ()=>assert.match(`${data}`, /37\.79264\d*/i))
		}
	})
	childProcess.stderr.on('data', data => {
		assert.fail(data)
	})
	
	await new Promise((r,x)=>{
		childProcess.on("close", ()=>setTimeout(r, 20))
	})
	console.log("end")
})
