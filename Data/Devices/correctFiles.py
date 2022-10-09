import sys

fileName = sys.argv[1]
fileName = fileName.split(".")[0]

search_text = ","

replace_text = " "

with open(fileName+".csv", 'r') as file:
	data = file.read()
	data = data.replace(search_text, replace_text)

with open(fileName+".dat", 'w') as file:
	file.write(data)

print("Done!")
