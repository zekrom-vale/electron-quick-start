// Tests the platform switch
const path = require('path')
process.env.NODE_ENV = process.platform
const config = require("config")
const assert = require('assert/strict');
const {it, describe, after} = require('node:test');

const execPath = path.normalize(
			path.join(__dirname, config.get("R.path.portable"))
	)
//Platform switch
const MACOS = "darwin"
const WINDOWS = "win32"
const LINUX = "linux"
const main = require("../main.js")

console.log("Platform Tests")
var home=path.join(__dirname, config.get("R.path.home"))


describe("Platform", ()=>{
	after(()=>{
		setTimeout(()=>{
			main.cleanUpApplication()
			app.exit(0)
		}, 20)
	})
	it( "MacOS", ()=>{
		assert.match(main.platform(MACOS, 1), /sed -i\s+""\s*'s!R_HOME_DIR=\.\*\$!R_HOME_DIR="[^!@#$%^&*()]+"\!'\s+[^!@#$%^&*()]+/)
	})
	it( "Windows", ()=>{
		assert.ifError(main.platform(WINDOWS))
	})
	it( "Linux", ()=>{
		assert.match(main.platform(LINUX, 1), /sed -i 's!R_HOME_DIR=\.\*\$!R_HOME_DIR="[^!@#$%^&*()]+"\!'\s+[^!@#$%^&*()]+/)
	})
	it( "SUN", ()=>{
		assert.throws(()=>main.platform("SUN"))
	})
})
