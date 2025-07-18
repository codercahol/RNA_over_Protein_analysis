module RnPAbxAnalysis

using DataFrames
using XLSX
using Plots
using Statistics
using LaTeXStrings
using LsqFit
using InvertedIndices
using Colors

function read_data(file)
    gc = DataFrame(XLSX.readtable(file, "Growth_od600"))
	protein = DataFrame(XLSX.readtable(file, "Protein_od555"))
	rna = DataFrame(XLSX.readtable(file, "RNA_od260"))
    bsa_standard = DataFrame(XLSX.readtable("BSA calibration curve.xlsx", "Sheet1"))
    
    xf = XLSX.readxlsx(file)
	sheets = XLSX.sheetnames(xf)
    if "RNA_OD600_wash_loss" in sheets
        rna_OD_wash_loss = DataFrame(XLSX.readtable(file, "RNA_OD600_wash_loss"))
    else
        rna_OD_wash_loss = DataFrame()
    end
	return gc, protein, rna, bsa_standard, rna_OD_wash_loss
end

function generate_growth_curve_plots(gc)
    # OD600 growth curve plot
    growth_plot = plot(gc.elapsed_mins, gc.OD600, label="")
	scatter!(
        gc.elapsed_mins, gc.OD600, label="", ylabel="OD600", 
        xlabel="Minutes Elapsed", title="Growth Curve"
    )

    # instantaneous growth rate plot
	plus_step = gc.OD600[2:end]
	minus_step = gc.OD600[1:end-1]
	t_plus = gc.elapsed_mins[2:end]
	t_minus = gc.elapsed_mins[1:end-1]
	dt = (t_plus - t_minus) / 60
	inst_GR = (plus_step - minus_step)./dt
	igr_plot = scatter(
        gc.elapsed_mins[2:end], inst_GR, label="", 
        ylabel="Instantaneous Growth Rate (1/hr)", xlabel="Elapsed Mins"
    )

    return growth_plot, igr_plot
end

function calculate_protein_OD555_conversion_to_biomass(bsa_standard)
    ODs::Vector{Float64} = bsa_standard[:,"OD (avg)"]
	prot_ladder::Vector{Float64} = bsa_standard[:,"conc. (ug/ml)"]

	# linear regression
	model(x, p) = p.*x
	X = ODs
	Y = prot_ladder
	fit = curve_fit(model, X, Y, [0.0])
	beta = fit.param[1]
    return beta
end

function raw_protein_mass(protein, beta)
    # 2x factor due to calibration curve of Plate reader to OD555 in quartz cuvette
	prot_effective_OD = skipmissing(protein.normed_PR) .* 2
	test_tube_vol = 1.5 # mL
	final_sample_vol = 0.18 # mL # vol of calibration curve ref.s
	prot_mass = prot_effective_OD .* beta * final_sample_vol
    return prot_mass
end

function protein_mass_fraction(prot_mass, gc, protein)
    biomass_per_OD600 = 450 # ug/mL
    test_tube_vol = 1.5 # mL
	elapsed_mins = collect(skipmissing(protein.elapsed_mins))
    protein_OD_loss = gc.Protein_supernatant_OD
	# in units OD*mL
    # biomass present at end of experiment
	biomass_for_prot = (gc.OD600 - protein_OD_loss) * test_tube_vol * biomass_per_OD600
	biomass_for_prot = filter(!ismissing, biomass_for_prot)
	spec_prot_conc = prot_mass ./ biomass_for_prot
	return spec_prot_conc
end

norm(x) = x[2:end] .- x[1]

function RNA_mass_fraction(gc, rna, rna_OD_wash_loss)
    rna_OD_loss = collect(skipmissing(gc.RNA_supernatant_OD))
    if !isempty(rna_OD_wash_loss)
        rna_secondary_OD_wash_loss = (
            norm(rna_OD_wash_loss.A) + norm(rna_OD_wash_loss.B) +
            norm(rna_OD_wash_loss.C) + norm(rna_OD_wash_loss.D)
        ) /4 * 1.2  # 1.2 mL is the total volume of washes, converted to extensive units (OD*mL)
        # 1.5 mL is the volume basis of original sample
        # converting back to intensive units
        rna_secondary_OD_wash_loss = rna_secondary_OD_wash_loss / 1.5

        rna_OD_loss = rna_OD_loss + rna_secondary_OD_wash_loss
    end

    biomass_per_OD600 = 450 # ug/mL
    normed_rna = collect(skipmissing(rna.OD260_normed))
    OD600 = gc.OD600[Not(ismissing.(gc.RNA_supernatant_OD))]
    biomass_for_rna = (OD600 - rna_OD_loss) * biomass_per_OD600 # mg/mL
	spec_rna = normed_rna * 31 ./ biomass_for_rna
    return spec_rna
end


const rRNA_per_totalRNA = 0.85
const RtoP = 1.0 / 1.7 # MDa protein / MDa RNA in a ribosome
const RNA_cal = rRNA_per_totalRNA * RtoP #* RNA_AA_mass_ratio

function phi_R(spec_rna, spec_prot_conc)
	R_prot = RNA_cal .* spec_rna
	phiR = R_prot ./ spec_prot_conc
    return phiR
end

function convert_RNA2Protein_ratio_to_phiR(ratio)
    """ data from Dai et el"""
    return RNA_cal * ratio
end

# is the total RNA per biomass consistent with the expected value from steady states at initial and final values?
# aim to have a comparison and a draft about this out really soon -- just want to have this done
# add one more small piece to manuscript about simulating fluctuating environments
# come with model comparison and overall structure of the manuscript ready to discuss on friday
# discuss manuscript overall structure on friday. Ship the manuscript by end of Jan

function RNA_mass_fraction(file::String; skip_invalid=false)
    gc, protein, rna, bsa_standard, rna_OD_wash_loss = read_data(file)
    growth_plot, igr_plot = generate_growth_curve_plots(gc)
    beta = calculate_protein_OD555_conversion_to_biomass(bsa_standard)
    prot_mass = raw_protein_mass(protein, beta)
    spec_prot_conc = protein_mass_fraction(prot_mass, gc, protein)
    spec_rna = RNA_mass_fraction(gc, rna, rna_OD_wash_loss)
    elapsed_mins = collect(skipmissing(protein.elapsed_mins))
    return spec_rna, elapsed_mins
end

function phi_R(file::String; skip_invalid=false)
    gc, protein, rna, bsa_standard, rna_OD_wash_loss = read_data(file)
    growth_plot, igr_plot = generate_growth_curve_plots(gc)
    beta = calculate_protein_OD555_conversion_to_biomass(bsa_standard)
    prot_mass = raw_protein_mass(protein, beta)
    spec_prot_conc = protein_mass_fraction(prot_mass, gc, protein)
    spec_rna = RNA_mass_fraction(gc, rna, rna_OD_wash_loss)
    elapsed_mins = collect(skipmissing(protein.elapsed_mins))
    if skip_invalid
        valid_idx = spec_prot_conc .+ spec_rna .<= 1 .&& spec_rna .>= 0.05
        valid_rna = spec_rna[valid_idx]
        valid_prot = spec_prot_conc[valid_idx]
        phiR = phi_R(valid_rna, valid_prot)
        ts = elapsed_mins[valid_idx]
    else
        phiR = phi_R(spec_rna, spec_prot_conc)
        ts = elapsed_mins
    end
    return phiR, ts
end

# module end
end