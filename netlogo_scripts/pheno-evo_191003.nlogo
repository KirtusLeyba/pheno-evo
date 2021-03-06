extensions [palette] ;; this allows us to use fancy color palettes

;;;;;;;;; variables

globals [
  ;; globals that won't change
  ;; we could even eventually delete the variables and write the values into the code...
  death-threshold ;; damage level below which a turtle dies ;; set to 0.0000001
  toxin-limit ;; concentration of toxin in a patch below which it goes to 0 ;; set to 0.0000001
  basal-growth-rate ;; probability of replicating at each timestep ;; set to 1
  tradeoff ;; constant by which to multiply growth rate as fitness tradeoff for having high degradation rate ;; set to 1
  alpha ;; defines shape of tradeoff of fitness and toxin degradation ;; set to 1
  pulsing? ;; boolean for pulsing ;; set to TRUE

  ;; globals we might end up using/changing
  switch-step ;; amount by which phenotype switching probabilty is allowed to mutate ;; for now, 0.1
  env-step ;; amount by which environmental response is allowed to mutate ;; for now, 0.1

  ;; globals that exist, but that we don't fiddle with
  ticks-to-pulse ;; number of ticks before next pulse of toxin; this acts as a counter

  ;; globals that are controlled by sliders/buttons... for now, anyway
  ; toxin-conc ;; initial concentration of toxin to add ;; set to 1
  ; n ;; number of turtles to initiate population
  ; initial-switch-rate ;; probability of switching phenotype (by any means) ;; may eventually be evolved
  ; initial-env-response ;; if switching, how important is the environmental signal? ;; may eventually be evolved
  ; diff-rate ;; amount of toxin that diffuses from one cell to the next at each timestep ;; important aspect of experiment
  ; env-noise ;; variance around the environmental signal that indicates how much toxin is in a patch ;; important aspect of experiment
  ; color-code ;; whether to color-code by growth rate or by phenotype. this is for model exploration only-- not for the long run
  ; pulse-rate ;; number of timesteps between antibiotic pulses ;; for now, 20... but should depend on the spatial structure of the pulse!
  ; mutation-rate ;; chance of mutating at each timestep

]


patches-own[
  toxin ;; toxin concentration on patch
  signal ;; how much toxin the patch *says* it has
]


turtles-own[
  mean1 ;; mean of 1st gamma distribution comprising genotype
  var1 ;; variance of 1st gamma distribution ;; if we're not evolving this, we could get rid of the variable
  wt1 ;; weight of 1st gamma distribution ;; note that the weight of the 2nd gamma distribution is 1 - wt1
  mean2 ;; mean of 2nd gamma distribution
  var2 ;; variance of 2nd gamma distribution ;; same as for var1
  switch-rate ;; frequency of switching phenotype
  degrade-rate ;; phenotypic rate at which cell degrades toxin, drawn from the genotype distribution
  degrade-energy ;; amount of energy being spent on degradation (updated whenever phenotype is updated)
  response-error ;; genotypically encoded noise in the individual's phenotypic response to the environment
  response-energy ;; amount of energy being spent on responding to the environment (updated whenever a response is made)
  growth-rate ;; this will be calculated from basal growth rate and degradation rate, accounting for fitness tradeoff
  health ;; health/damage level
  barcode ;; unique identifier of the cell lineage
]


;;;;;;;;; small functions

to poison ;; initiate toxin in one area
   ask patches[
  ; ask patches [if pxcor mod 2 = pycor mod 2 [
  ; ask patches with [(pxcor > 20) and (pxcor < 30)] [ ;; swath of poison down the middle
  ; ask patch random max-pxcor random max-pycor [ ;; choose random spot
    set toxin toxin + toxin-conc ;; add the "toxin-conc" amount to that patch
    ]
  ; ] ;; note: if this isn't included in the "go" function, need to make sure to update patch color and signal afterward.
end


to dilute ;; simulate transfer to fresh medium
  if dilute-rate != 1 [
    let newpop-num round (count turtles / dilute-rate) ;; figure out how many turtles left after dilution
    let turtles-to-save n-of newpop-num turtles ;; choose that many random turtles from the population
    ask turtles [ ;; ask turtles whether they're part of the population to be saved
      if not member? self turtles-to-save
      [die] ;; if they aren't, ask them to die
    ]
    ask turtles [
      if count patches with [ not any? turtles-here] > 0 [
        move-to one-of patches with [ not any? turtles-here ] ;; redistribute the survivors randomly across the environment
      ]
    ]
  ]
end


to color-turtles ;; color-code cells either by phenotype or by growth rate, as desired
  if color-code = "growth-rate" [
    ask turtles[
      set color palette:scale-gradient [[197 27 125][233 163 201][161 215 106][77 146 33]] growth-rate 0 1
    ]
  ]
  if color-code = "degrade-rate" [
    ask turtles[
      set color palette:scale-gradient [[165 0 38][255 255 191][49 54 149]] degrade-rate 0 1
    ]
  ]
  if color-code = "barcode" [
    ask turtles[
      set color palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 11 barcode 0 (n * 10)
    ]
  ]
end


;;;;;; important functions

to make-environment
  set death-threshold 0.0000001
  set toxin-limit 0.0000001
  set basal-growth-rate 1
  set tradeoff 1
  set alpha 1
  set pulsing? TRUE
  set switch-step 0.1
  set env-step 0.1
  ask patches[ ;; initiate patches
    set toxin 0
    set signal toxin
    set pcolor white
  ]
end


to add-cells
  crt n [ ;; make cells
    setxy random max-pxcor random max-pycor ;; scatter randomly
    set size 1
    set shape "square" ;; this is unnecessary but helps identify the ancestral cells in case that's interesting
    set switch-rate initial-switch-rate
    set response-error initial-response-error
    set degrade-rate [signal] of patch-here + random-float response-error - 0.5 * response-error
    if degrade-rate < 0 [set degrade-rate 0]
    if degrade-rate > 1 [set degrade-rate 1]
    set response-energy (1 - abs ([signal] of patch-here - degrade-rate)) / 2 ;; not sure this is the best way to calculate this...
    set degrade-energy degrade-rate / 2
    set growth-rate basal-growth-rate - degrade-energy - response-energy
    set health 1 ;; everyone starts healthy
    set barcode random n * 10 ;; assign a random barcode to help identify the lineage
  ]
    color-turtles ;; color-code turtles by growth rate or phenotype
end


to setup
  clear-all
  make-environment
  add-cells
  reset-ticks
end


;;;;;;; making the model go

to go
  ;;;; turtle stuff
  ask turtles[

    ;;; switch phenotype. do this first so that other behaviors can work according to phenotype
    if random-float 1 < switch-rate[ ;; first, decide whether to switch. if switching,
      ;; set degrade-rate to match environmental signal, plus a stochastic amount of wrongness determined by genetic propensity for wrongness
      set degrade-rate [signal] of patch-here + random-float response-error - 0.5 * response-error
      ifelse degrade-rate < 0 [set degrade-rate 0]
      [if degrade-rate > 1 [set degrade-rate 1]]
      set response-energy (1 - abs ([signal] of patch-here - degrade-rate)) / 2 ;; the energy spent on sensing can be calculated from signal and degrade-rate
      set degrade-energy degrade-rate / 2
    ]

    ;;; effect of toxin on turtles and vice versa
    ;;; here, toxin poisons turtles before turtles get a chance to degrade. we could change that if we want.
    ifelse [toxin] of patch-here > 0 [ ;; check whether toxin is present
      set health (health - (toxin)) ;; subtract health by amount equal to toxin level. note there's no coefficient on toxin
      ifelse health < death-threshold [ ;; then, check whether health is below threshold
        die ;; if it is, die
      ]
      [let dr degrade-rate ;; if health > death-threshold, start degrading toxin
        ask patch-here [
          set toxin (toxin - dr) ;; degrade some of it by according to phenotype. note there's no coefficient on dr
          if toxin < toxin-limit [
            set toxin 0
          ]
        ]
      ]
    ]
    [set health 1];; if [toxin] of patch-here = 0, recover health completely

    ;;; reproduction and mutation
    ;;; first, figure out growth rate based on how much energy was spent on other activities
    set growth-rate basal-growth-rate - degrade-energy - response-energy
    ;;; then, make new cells
    if (count neighbors with [not any? turtles-here]) > 0 [ ;; only reproduce if there's space
      if random-float 1 < growth-rate [ ;; higher growth rate -> higher probability of reproducing
        ;; store the values that are heritable - for reproduction
        let new-switch-rate [switch-rate] of self
        let new-response-error [response-error] of self
        let new-degrade-rate [degrade-rate] of self
        let new-barcode [barcode] of self
        let new-growth-rate [growth-rate] of self
        ;; mutate, if necessary
        if mutation-rate != 0 [ ;; only bother with this if mutation rate isn't 0
          if random-float 1 < mutation-rate [ ; do a random draw to see whether it's time to mutate
            set new-switch-rate [switch-rate] of self + random-float switch-step - 0.5 * switch-step
            set new-response-error [response-error] of self + random-float env-step - 0.5 * env-step
          ]
        ]
        ;; now reproduce
        let growth-space patch-set neighbors with [not any? turtles-here] ;; find patches with no turtles
        ask one-of growth-space [ ;; choose one of the empty patches to sprout
          sprout 1 [ ;; below are all the properties of the newly sprouted cell
            set shape "bacterium"
            set switch-rate new-switch-rate
            set response-error new-response-error
            set degrade-rate new-degrade-rate
            set barcode new-barcode
            set growth-rate new-growth-rate
            set health 1
          ]
        ]
      ]
    ]
  ]

  ;;;; dilute the population whenever the universe fills up
  if count turtles >= count patches [
    dilute
  ]

  ;;;; to make toxin get added in periodic pulses
  ;if pulsing? [ ;; for now, got rid of the "if" statement because we're always pulsing
    if ticks-to-pulse <= 0 [ ;; check the ticks-to-pulse counter. if it's at zero...
      poison ;; add a toxin pulse...
      set ticks-to-pulse pulse-rate ;; and reset the ticks-to-pulse counter
    ]
    set ticks-to-pulse (ticks-to-pulse - 1) ;; if the counter isn't yet at zero, subtract one
  ;]

  ;;;; patches: manage toxin levels
  diffuse toxin diff-rate ;; diffuse toxin
  ask patches[
    if toxin < toxin-limit [set toxin 0] ;; if toxin goes below threshold, make it zero
    ;; signal to the world how much toxin each patch has (with some noise):
    set signal toxin + random-float env-noise - 0.5 * env-noise
    ifelse signal < 0 [set signal 0]
    [if signal > 1 [set signal 1]]
    ;; color-code patches by toxin concentration: grayscale on log scale with low = light and high = dark:
    ifelse toxin > 0 [set pcolor scale-color gray (log toxin 10) 1 (log toxin-limit 10)]
    [set pcolor white]
  ]

  ;;;; color-code turtles by phenotype or growth rate
  color-turtles

  ;;;; time elapses
  tick

end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
728
529
-1
-1
10.0
1
10
1
1
1
0
0
0
1
0
50
0
50
1
1
1
ticks
30.0

BUTTON
13
10
79
43
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
12
493
75
526
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
13
44
63
104
n
200.0
1
0
Number

PLOT
730
11
890
131
total turtle population
time
turtles
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles"

PLOT
892
10
1052
130
average toxin per patch
time
toxin
0.0
100.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [toxin] of patches"

PLOT
732
374
892
494
degradation rate phenotype
degradation rate
turtles
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" "histogram [degrade-rate] of turtles"

SLIDER
13
107
185
140
initial-switch-rate
initial-switch-rate
0
1
0.0
.1
1
NIL
HORIZONTAL

SLIDER
12
144
198
177
initial-response-error
initial-response-error
0
1
1.0
.1
1
NIL
HORIZONTAL

SLIDER
12
217
184
250
diff-rate
diff-rate
0
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
12
181
184
214
env-noise
env-noise
0
1
1.0
0.1
1
NIL
HORIZONTAL

PLOT
891
132
1051
252
patch signal v. toxin
toxin
signal
0.0
1.0
0.0
1.0
true
false
"" "clear-plot"
PENS
"default" 0.01 2 -16777216 true "" "ask patches [plotxy toxin signal]"

CHOOSER
12
254
150
299
color-code
color-code
"growth-rate" "degrade-rate" "barcode"
1

PLOT
731
253
891
373
turtle responses to patches
signal
deg-rate
0.0
1.0
0.0
1.0
true
false
"" "clear-plot"
PENS
"default" 0.01 2 -16777216 true "" "ask turtles [plotxy [signal] of patch-here degrade-rate]"

PLOT
893
375
1053
495
turtle growth rates
growth rate
turtles
0.0
1.1
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" "histogram [growth-rate] of turtles"

PLOT
894
254
1054
374
turtle energy expenses
degrade
response
0.0
0.4
0.0
0.4
true
false
"" "clear-plot"
PENS
"default" 0.01 2 -16777216 true "" "ask turtles [plotxy degrade-energy response-energy]"

PLOT
730
132
890
252
toxin in patches
toxin
patches
0.0
1.0
0.0
500.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" "histogram [toxin] of patches"

INPUTBOX
64
44
128
104
pulse-rate
20.0
1
0
Number

SLIDER
12
303
184
336
mutation-rate
mutation-rate
0
1
1.0
0.1
1
NIL
HORIZONTAL

PLOT
1055
376
1215
496
genotypic switch rates
NIL
NIL
0.0
1.1
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" "histogram [switch-rate] of turtles"

PLOT
1056
255
1216
375
genotypic response error
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" "histogram [response-error] of turtles"

INPUTBOX
130
45
197
105
toxin-conc
0.5
1
0
Number

SLIDER
10
345
182
378
dilute-rate
dilute-rate
1
100
100.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

bacterium
true
0
Polygon -7500403 true true 135 30 165 30 195 45 195 195 165 210 135 210 105 195 105 45
Line -7500403 true 150 210 150 210
Line -7500403 true 150 210 150 210
Line -7500403 true 150 210 150 240
Line -7500403 true 150 240 135 255
Line -7500403 true 135 255 150 285

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment_190830" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="10000"/>
    <metric>count turtles</metric>
    <metric>mean [toxin] of patches</metric>
    <metric>[toxin] of patches</metric>
    <metric>[degrade-rate] of turtles</metric>
    <enumeratedValueSet variable="toxin-conc">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-response-error">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-noise">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pulse-rate">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-switch-rate">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-code">
      <value value="&quot;degrade-rate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="diff-rate">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
