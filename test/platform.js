// Tests the platform switch
const {app, BrowserWindow} = require('electron')
const path = require('path')
process.env.NODE_ENV = process.platform
const config = require("config")
const assert = require('assert/strict');
const {it, describe, after} = require('node:test');

//Platform switch
const MACOS = "darwin"
const WINDOWS = "win32"
const LINUX = "linux"
const main = require("../main.js")

describe("Platform", ()=>{
	after(()=>{
		setTimeout(()=>{
			main.cleanUpApplication()
			app.exit(0)
		}, 20)
	})
	it( WINDOWS, ()=>{
		assert.ifError(main.platform(WINDOWS))
	})
	it( LINUX, ()=>{
		assert.match(main.platform(LINUX, 1), /sed -i 's!R_HOME_DIR=\.\*\$!R_HOME_DIR="[^!@#$%^&*()]+"!' [^!@#$%^&*()]+$/)
	})
	it( MACOS, ()=>{
		assert.match(main.platform(MACOS, 1), /sed -i\s+""\s*'s!R_HOME_DIR=\.\*\$!R_HOME_DIR="[^!@#$%^&*()]+"!' [^!@#$%^&*()]+$/)
	})
	it( "SUN (Other)", ()=>{
		assert.throws(()=>main.platform("SUN"))
	})
})
