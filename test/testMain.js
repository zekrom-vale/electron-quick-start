const assert = require('assert/strict');
const test = require('node:test');

console.log("Main Test")
test("Main Test", async (contex) =>{
	const main = require("../main.js")
	await new Promise((r,x)=>{setTimeout(r, 1*6000)})
	main.cleanUpApplication()
	console.log("end")
})