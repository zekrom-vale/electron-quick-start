module.exports = {
	start: async function (platform, arch, ignore){
		const path = require('path')
		const {execSync, exec} = require('child_process')
		const fs = require('fs-extra');
		
		process.env.NODE_ENV = platform
		const config = require("config")
		
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
	}
}