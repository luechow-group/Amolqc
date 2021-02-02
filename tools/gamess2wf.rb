#!/usr/bin/env ruby

=begin
 Copyright (C) 2013 Alexander Sturm
 Copyright (C) 2014, 2017 Kaveh Haghighi Mood

 SPDX-License-Identifier: GPL-3.0-or-later
=end

	def sgnblk(a)
		num=1
		for i in 1..a.size-1
					num+=1 if a[i-1]*a[i]<0
    end
    return num
	end

if ARGV.size < 2 then
	puts <<EOD
		Usage: gamess2wf.rb outfile punchfile [basis] [determinant cutoff]
		Note: this script only supports restricted HF/MCSCF calculations so far.
EOD
	exit
end
cutoff = 0
basis = "general"
if ARGV.size == 4  then
	basis  = ARGV[2].to_s
	cutoff = ARGV[3].to_f
	ARGV.pop(2)
elsif ARGV.size ==3
		basis  = ARGV[2].to_s
		ARGV.pop
end
charge = 0
mult = 1
ncore = 0
nelec = 0
mcscf = false
maxorb = nil
dets = []
csfdets = []
csfs = []
csf = []
natorbs = false
mos = []
geom = []

while ARGF.gets do
	ncore = $_.split[-1].to_i if $_.include? "NUMBER OF CORE ORBITALS"
	if ncore == 0
		then
		 ncore = $_.split.at(1).to_i if $_.include? "NFZC"
		 ncore = $_.split.at(1).to_i if $_.include? "NMCC"
		 corarr=Array(1..ncore)
	end

	if $_.include? "NUMBER OF ELECTRONS" then
		nelec = $_.split[-1].to_i
		maxorb ||= (nelec + 1) / 2
	end
	charge = $_.split[-1].to_i if $_.include? "CHARGE OF MOLECULE"
	mult = $_.split[-1].to_i if $_.include? "SPIN MULTIPLICITY"
	mcscf = true if $_.include? "SCFTYP=MCSCF"
	natorbs = true if $_.include? "FINCI=NOS"
	(mcscf = true; maxorb = 10000) if $_.include? "SCFTYP=NONE"

	if mcscf then
		start = $_[/^\s+STATE\s+1\s+ENERGY/]
		if start ... (~/^$/ or $_.include?("DONE WITH")) then
			if start then
				dets.clear
				(1..4).each do ARGF.gets end
			end
			alpha, _, beta, _, coeff, rest = $_.split
			next unless alpha and beta and coeff
			next if rest
			next if coeff.to_f.abs < cutoff
			det = [[], []]
			[alpha, beta].each_with_index do |conf, ci|
				det[ci].concat [*1..ncore]
				conf.chars.each_with_index do |char, idx|
					orb = idx + 1 + ncore
					det[ci] << orb if char == '1'
					maxorb = [maxorb, orb].max
				end
			end
			dets << [coeff.to_f, det.join(" ")]
		end

		if $_.include?("DETERMINANT CONTRIBUTION TO CSF") .. $_.include?("TOTAL NUMBER OF INTEGRALS") then
			if $_.include? " C(" then
				coeff = $_.split("=")[-1].to_f
				det = $_.split(":")[-1]
				alpha = []
				beta = []
				detint=[]
				det=det.gsub("-"," -")
				detint=det.split.reverse
				detint=detint.map{|elem| elem.to_i}

              if sgnblk(detint)== 2  then
                       if (detint[0]*detint[1]<0)
                        if detint[1]>detint[0]
                          coeff=-coeff
                          tmp=detint[1]
                          detint[1]=detint[0]
                          detint[0]=tmp
                        end
                      end
               else
                  while sgnblk(detint)>2
                  	for i in 1..detint.size-1
                  		if (detint[i-1]*detint[i]<0)
                  			if detint[i]>detint[i-1]
                  				coeff=-coeff
                  				tmp=detint[i]
                  				detint[i]=detint[i-1]
                  				detint[i-1]=tmp
                  			end
                  		end
                    end
                  end
               end
				detint.each do |el|
					el = el.to_i
					if el < 0 then
						el = -el
						beta << el
					else
						alpha << el
					end
					maxorb = [maxorb, el].max
				end
  				#beta.each{|element| print element," "}
					#print "\n"
        while alpha != alpha.sort
        	for i in 1..alpha.size-1
        		if (alpha[i]< alpha[i-1]) then
        			coeff=-coeff
        			tmp=alpha[i]
        			alpha[i]=alpha[i-1]
        			alpha[i-1]=tmp
        		end
        	end
        end
        while beta != beta.sort
        	for i in 1..beta.size-1
        		if (beta[i]< beta[i-1]) then
        			coeff=-coeff
        			tmp=beta[i]
        			beta[i]=beta[i-1]
        			beta[i-1]=tmp
        		end
        	end
        end
				if ncore !=0 then
					alpha.insert(0,*corarr)
					 beta.insert(0,*corarr)
				end
				csf << "  #{coeff} #{alpha.join(" ")} #{beta.join(" ")}"
			else
				csfdets << csf unless csf.empty?
				csf = []
			end
		end

		start = $_[/CSF\s+COEF/]
		if start .. ($_.include?("END OF CI-MATRIX DIAGONALIZATION") or ~/^$/) then
			(csfs.clear; next) if start
			words = $_.split
			csf = words[0].to_i
			coef = words[1].to_f
			next unless csf > 0 and coeff.abs > cutoff
			csfs << [coef, csfdets[csf - 1]]
		end
	end

	start = $_.start_with?(" $VEC") && (not natorbs or mos.empty?)
	stop = $_.include? "$END"
	if start ... stop then
		(mos.clear; next) if start
		mos << $_ unless stop
	end

	start = $_.include? "COORDINATES"
	if start ... ~/^$/ then
		next if start or $_.include? "CHARGE" or ~/^$/
		atom, _, x, y, z = $_.split
		geom << " #{atom} #{x} #{y} #{z}"
	end
end
puts "$general", "title='No title'", "charge=#{charge}, spin=#{mult}, geom=bohr",
     "evfmt=gms,basis=#{basis} , jastrow=none", "$end"
puts "$geom", "#{geom.size}", geom.map{|l|l.strip.capitalize}, "$end"

orbs = mos.select{|line|line.to_i <= maxorb}
puts "$mos", "#{orbs[-1].split[0]}", "", orbs, "$end"
if not dets.empty? and not csfs.empty? then
	$stderr.puts "WARNING! Found both CSFs and determinants!"
elsif dets.empty? and csfs.empty? then
	#$stderr.puts "WARNING! Found no CSFs/determinants!"
	#$stderr.puts "WARNING! Please Check CSFs/determinants!"
	puts "$csfs"
  puts " single restricted"
  puts "$end"
end

if not dets.empty?
dets = dets.sort_by{|c,d|-c.abs}
	puts "$dets", "#{dets.size}"
	0.upto(dets.length-1)  { |i| print dets[i][0]," ",dets[i][1],"\n"}
	puts	"$end"
end
if not csfs.empty?
	puts "$csfs", "#{csfs.size}"
	csfs.sort_by{|coef,_| -coef.abs }.each do |coef, dets|
		puts "#{coef} #{dets.size}"
		puts dets
	end
	puts "$end"
end

