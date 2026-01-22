#
# viral_surveillance
#
viral_surveillance:
	womtool validate --inputs viral_surveillance.json viral_surveillance.wdl
	miniwdl check viral_surveillance.wdl
	#womtool inputs viral_surveillance.wdl|jq -S -M --indent 4 > viral_surveillance.json

viral_surveillance_docu:
	wdl-aid viral_surveillance.wdl -o viral_surveillance.md
	womtool graph viral_surveillance.wdl > viral_surveillance.dot
	dot -Tpdf -o viral_surveillance.pdf viral_surveillance.dot
	dot -Tjpeg -o viral_surveillance.jpeg viral_surveillance.dot
	rm viral_surveillance.dot

run_viral_surveillance:
	miniwdl run --debug --dir test-viral_surveillance --cfg miniwdl_production.cfg --input viral_surveillance.json viral_surveillance.wdl

