const path = require('path')
process.env.NODE_ENV = process.platform
const config = require("config")
const assert = require('assert/strict');
const test = require('node:test');

const execPath = path.normalize(
			path.join(__dirname, config.get("R.path.portable"))
	)
//Platform switch
const MACOS = "darwin"
const WINDOWS = "win32"
const LINUX = "linux"
function platform(p){
	// let and const here will be discarded after this block, var will be kept
	let _i=''
	switch(p){
		case WINDOWS:
			return null
			break
		case MACOS:
			// Adapt to support MacOS sed
			_i=' ""'
			//Fall through as mac is linux like
		case LINUX:
			// Fix issue with R Home path by overriding the R sh script with the correct value
			let home=path.join(__dirname, config.get("R.path.home"))
			// Must use ! as / is an issue with paths
			return(`sed -i${_i} 's!R_HOME_DIR=.*$!R_HOME_DIR="${home}"!' ${execPath}`)
			break
		default:
			throw new Error(`Platform Error: Not on windows, linux, or macos. Got ${p}`)
			
	}
}
console.log("Platform Tests")
var home=path.join(__dirname, config.get("R.path.home"))
test( "Platform Test - MacOS", contex =>{
	assert.strictEqual(platform(MACOS), `sed -i "" 's!R_HOME_DIR=.*$!R_HOME_DIR="${home}"!' ${execPath}`)
})
test( "Platform Test - Windows", contex =>{
	assert.ifError(platform(WINDOWS))
})
test( "Platform Test - Linux", contex =>{
	assert.strictEqual(platform(LINUX), `sed -i 's!R_HOME_DIR=.*$!R_HOME_DIR="${home}"!' ${execPath}`)
})
test( "Platform Test - SUN", contex =>{
	assert.throws(()=>platform("SUN"))
})
console.log("end")