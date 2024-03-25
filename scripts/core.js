module.exports = {
	start: function (platform, arch, ignore = undefined){
		const path = require('path')
		const {execSync} = require('child_process')
		const fs = require('fs-extra');
		
		process.env.NODE_ENV = platform
		const config = require("config")
		if(ignore === undefined)ignore = config.get("app.ignore")
		if(!config.get("R.path.isPortable")){
			if(ignore=="") ignore = "R-Portable-*"
			else ignore = `R-Portable-*|${ignore}`
		}
		
		var name=config.get("app.name")
		var out=config.get("app.out")
		var buildPath=path.join(process.cwd(), out, `${name}-${platform}-${arch}`)
		console.log(buildPath)

		var run=`electron-packager .
${name}
--overwrite
--platform=${platform}
--arch=${arch}
--out=${out}
--icon=${config.get("app.icon")}
--prune=true
--version-string.CompanyName=${config.get("app.CompanyName")}
--version-string.FileDescription=${config.get("app.FileDescription")}
--version-string.ProductName=${config.get("app.ProductName")}
--ignore="${ignore}|\\.git.*|\\.Rhistory|.*\\.Rproj|scripts|${out}"`
		run=run.replace(/[\n\r]+/g, " ")
		console.log(run)
		execSync(run)

		fs.copySync(path.join(buildPath, "resources", "app", "config"), path.join(buildPath, "config"))
	},
	
	test: async function (platform, arch){
		const path = require('path')
		const {execSync} = require('child_process')
		const fs = require('fs-extra');
		
		process.env.NODE_ENV = platform
		const config = require("config")

		var name=config.get("app.name")
		var out=config.get("app.out")
		var buildPath=path.join(process.cwd(), out, `${name}-${platform}-${arch}`)
		console.log(buildPath)
		
		execSync(`timeout 30 ${path.join(buildPath, name)}`)
	}
}