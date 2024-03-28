const assert = require('assert/strict');
const test = require('node:test');

console.log("Main Test")
test("Main Test", async (contex) =>{
	const main = require("../main.js")
	await new Promise((r,x)=>{setTimeout(r, 6*1000)})
	main.cleanUpApplication()
	console.log("end")
	assert.ok(true)
})