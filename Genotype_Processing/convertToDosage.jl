#Takes hapmap file line-by-line and outputs csv of allele dosage coded by major/minor allele frequency

using StatsBase

inFileName = ARGS[1]

function convert_SNP(callArray)
	combAlleleStr = string(filter(!=("NN"), callArray)...)	

	alleleCounts = countmap(combAlleleStr)

	majAllele = findmax(alleleCounts)[2]
	minAllele = findmin(alleleCounts)[2]

	callArray[callArray.==string(majAllele, majAllele)] .= "0,"
	callArray[callArray.==string(minAllele, minAllele)] .= "2,"
	callArray[callArray.==string(majAllele, minAllele)] .= "1,"
	callArray[callArray.==string(minAllele, majAllele)] .= "1,"
	callArray[callArray.=="NN"] .= "NA,"

	return(callArray)
end

print(string("Processing file ", inFileName, "\n"))

#Create file to write finished data to.

outFileName = replace(inFileName, ".hmp.txt" => "_dosage.csv")
outFile = open(outFileName, "w")

open(inFileName) do inFile
	for SNPLine in eachline(inFile)
		lineArray = split(SNPLine, "\t")
		SNPName = lineArray[1]

		#Delete our filler lines
		deleteat!(lineArray, 1:24) 

		#Is this the header line?
		if SNPName == "\"rs#\""
			newLine = string(lineArray...)
			newLine = replace(newLine, "\""=>"")
			newLine = replace(newLine, "UX"=>",UX")
			write(outFile, string("SNP", newLine, "\n"))
		else
			callCounts = countmap(lineArray)

			#Remove NAs from counts
			delete!(callCounts, "NN")

			if length(callCounts) <= 3
				dosageArray = convert_SNP(lineArray)
				newLine = chop(string(SNPName, ",", dosageArray...))
				write(outFile, string(newLine, "\n"))
			else
				print(string("SNP ", SNPName, " discarded with ",
					     length(callCounts), " classes."))
			end
		end
	end
end


close(outFile)
