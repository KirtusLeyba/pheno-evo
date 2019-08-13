extensions [palette] ;; this allows us to use fancy color palettes

globals [
  ;; patch globals
  max-density ;; maximum number of cells allowed per patch
  initial-food ;; initial amount of food per patch
  infinite-food? ;; should food be always replenished? (1,0)
  diff-rate ;; amount food and toxin that diffuses from one cell to the next at each timestep
  toxin-limit ;; lower limit of toxin below which it goes to 0
  ;; turtle globals
  growth-threshold ;; damage level below which a turtle can't replicate
  death-threshold ;; damage level below which a turtle dies
  repair-rate ;; increase in health per level time if there's no toxin
  pulsing? ;; boolean for pulsing
  ; turbidostat? ;; boolean for diluting
  ticks-to-pulse ;; number of ticks before next pulse of toxin
  turbidostat-ticks ;; number of ticks after diluting before pulsing
  ;; globals that are set by buttons
  ; n ;; number of turtles to initiate population
  ; basal-growth-rate ;; probability of replicating at each timestep (per timestep)
  ; initial-mu ;; initial value that makes up list of mu-values
  initial-sig ;; initial value that makes up list of sig-values
  initial-switch-prob ;; initial probability of switching phenotype
  alpha ;; exponent for relationship between growth rate and toxin degradation
  ; tradeoff ;; constant by which to multipleate growth rate as fitness tradeoff for having high degradation rate
  ; mu-step ;; amount by which mu is allowed to mutate each step
  sig-step ;; amount by which sigma is allowed to mutate each step
  ; switch-step;; amount by which phenotype switching probabilty is allowed to mutate
  ; mutation-rate ;; chance of mutating at each timestep
  ; damage-rate ;; how much damage a cell accrues due to toxin  (percent per timestep)
  heritable? ;; whether daughter cells inherit parental phenotype directly (1,0)
  ; toxin-conc ;; initial concentration of toxin to add
  ; dilute-rate ;; factor by which population is diluted when necessary
]


;; A note on the genotype:
;; it's the sum of 5 equally weighted uniform distributions,
;; each with mean mu and width sigma
;; so it extends from mu - 0.5 * sigma to mu + 0.5 * sigma
;; and must all fall between 0 and 1.
;; To set the phenotype,
;; we draw once from each of the distributions,
;; and then normalize by dividing by 5.
;; If later we decide to weight the distributions differently,
;; we will have to change the probability with which we draw from each distribution.


patches-own[
  toxin ;; toxin concentration on patch
  food ;; food concentration on patch
]


turtles-own[
  mu-vals ;; means of uniform distributions comprising genotype (a list of length 5)
  sig-vals ;; widths of uniform distributions comprising genotype (a list of length 5)
  switch-prob ;; probability of switching phenotype
  degrade-rate ;; rate at which cell degrades toxin, drawn from the genotype distribution
  health ;; health/damage level
  order-in-stack ;; where it is in the stack of cells on the patch (if multiple cells per patch are allowed)
  barcode ;; unique identifier of the cell lineage
  growth-rate ;; this will be calculated from basal growth rate and degradation rate, accounting for fitness tradeoff
]


to color-turtles ;; this determines the color scale for turtles, according to phenotype: low degradation ability = red, medium = yellow, high = blue
  set color palette:scale-gradient [[165 0 38][255 255 191][49 54 149]] degrade-rate 0 1
end


to color-patches ;; this determines the color scale for toxin concentration: grayscale on log scale with low = light and high = dark
  ifelse toxin > 0 [set pcolor scale-color gray (log toxin 10) 1 (log toxin-limit 10)]
  [set pcolor white]
end


to make-environment
  set max-density 1
  set initial-food 1
  set infinite-food? 1
  set diff-rate 0.2
  set toxin-limit 0.0000001
  set repair-rate 0.01
  set growth-threshold 0.5
  set death-threshold 0.1
  set pulsing? false ;; this is the default unless button says otherwise
  set heritable? true
  ask patches[ ;; initiate patches
    set toxin 0
    set food initial-food
    set pcolor white
  ]
end

to environment-1
  ;; Daniel's parameters
  ;; override variables determined by buttons; use these presets
  set basal-growth-rate 0.01
  set damage-rate 0.38
  set heritable? true
  set dilute-rate 5
  set toxin-conc 0.06
  set alpha 1
  set tradeoff 0.43
  set mutation-rate 0.005
  set mu-step 1
  set sig-step 0.25
  set initial-switch-prob 0.005
  set initial-sig 0.28
  set initial-mu 0.43
  set switch-step 0.1
  set pulse-rate 10
  set n 200
  set turbidostat? true
end

to environment-2
  ;; Kirtus's parameters
  ;; override variables determined by buttons; use these presets
  set basal-growth-rate 0.879
  set damage-rate 0.7
  set heritable? true
  set dilute-rate 5
  set toxin-conc 0.6
  set alpha 1
  set tradeoff 0.71
  set mutation-rate 0.512
  set mu-step 0.51
  set sig-step 0.2
  set initial-switch-prob 0.005
  set initial-sig 0.3
  set initial-mu 0.15
  set switch-step 0.2739
  set pulse-rate 53
  set n 20
  set turbidostat? false
  set pulsing? false
end

to environment-3
  ;; Jess's parameters
  ;; override variables determined by buttons; use these presets
  set basal-growth-rate 0.5
  set damage-rate 1
  set heritable? false
  set dilute-rate 10
  set toxin-conc 0.1
  set alpha 1
  set tradeoff 0.9
  set mutation-rate 0.005
  set mu-step 1
  set sig-step 0.25
  set initial-switch-prob 0.5
  set initial-sig 0.28
  set initial-mu 0.43
  set switch-step 1
  set pulse-rate 10
  set n 200
  set turbidostat? true
end

to add-cells
  crt n [ ;; make cells
    setxy random max-pxcor random max-pycor ;; scatter them randomly
    set size 1
    set shape "square" ;; this is unnecessary but helps identify the ancestral cells in case that's interesting
    if initial-mu + 0.5 * initial-sig > 1 [set initial-sig 2 * (1 - initial-mu)] ;; this corrects sigma in case the values set on the sliders don't adhere to our requirement
    if initial-mu - 0.5 * initial-sig < 0 [set initial-sig 2 * initial-mu] ;; that the distribution centered on mu with width sigma be between 0 and 1
    ;set mu-vals [0.1 0.3 0.5 0.7 0.9]
    set mu-vals n-values 5 [initial-mu] ;; the population is initialized with a set of 5 identical mu and sigma values, set by the sliders
    set sig-vals n-values 5 [initial-sig]
    foreach [0 1 2 3 4] [ [x] -> ; loop through all mu-values and all sig-values. at each list index between 0 and 4...
      if item x mu-vals > (1 - 0.5 * item x sig-vals) [
        set mu-vals replace-item x mu-vals (1 - 0.5 * item x sig-vals)]
      if item x mu-vals < 0.5 * item x sig-vals[
        set mu-vals replace-item x mu-vals (0.5 * item x sig-vals)]
    ]
    set switch-prob initial-switch-prob ;; set initial probability of switching phenotype
    set degrade-rate ( ;; determine phenotype (toxin degration capacity) by drawing from the 5-part genotypic distribution
      0.2 * (random-float item 0 sig-vals + (item 0 mu-vals - 0.5 * item 0 sig-vals) +
      random-float item 1 sig-vals + (item 0 mu-vals - 0.5 * item 1 sig-vals) +
      random-float item 2 sig-vals + (item 0 mu-vals - 0.5 * item 2 sig-vals) +
      random-float item 3 sig-vals + (item 0 mu-vals - 0.5 * item 3 sig-vals) +
      random-float item 4 sig-vals + (item 0 mu-vals - 0.5 * item 4 sig-vals)))
    color-turtles
    set health 1 ;; everyone starts healthy
    set order-in-stack 1 ;; everyone starts at the ground level
    set barcode random n * 10 ;; assign a random barcode to help identify the lineage
    set growth-rate basal-growth-rate * (1 - tradeoff * degrade-rate ^ alpha) ;; determine this cell's growth rate based on its degradation phenotype
  ]
end



to setup
  clear-all
  make-environment
  add-cells
  reset-ticks
end


to poison ;; initiate toxin in one spot - this is mostly just to test diffusion
  ask patch random max-pxcor random max-pycor [ ;; choose random spot
    set toxin toxin + toxin-conc ;; add the "toxin-conc" amount to that patch
    color-patches
  ]
end



to pulse ;; add toxin in a pulse across entire environment
  ask patches[
    set toxin toxin + toxin-conc ;; add the "toxin-conc" amount to that patch
    color-patches
  ]
end



to dilute ;; simulate transfer to fresh medium
  ; export-world (word "results " date-and-time ".csv") ;; note that this results in a rather large file - commenting it out for now
  let newpop-num round (count turtles / dilute-rate) ;; figure out how many turtles left after dilution
  let turtles-to-save n-of newpop-num turtles ;; choose that many random turtles from the population
  ask turtles [ ;; ask turtles whether they're part of the population to be saved
    if not member? self turtles-to-save
    [die] ;; if they aren't, ask them to die
  ]
  ask turtles [
    setxy random max-pxcor random max-pycor ;; redistribute the survivors randomly across the environment
  ]
  ; make-environment ;; reset environmental conditions (toxin, food, etc)
end


to start-pulse-rate ;; this allows a button to set off a periodic pulse of antibiotics
  set pulsing? true ;; once pulsing? is true, pulses will be part of the "go" command
  set ticks-to-pulse pulse-rate ;; counter for determining when to pulse next
end


to stop-pulse-rate ;; this allows a button to set off a periodic pulse of antibiotics
  set pulsing? false ;; once pulsing? is true, pulses will be part of the "go" command  set ticks-to-pulse pulse-rate ;; counter for determining when to pulse next
end



to go
  ; if not any? turtles [stop] ;; once the population is extinct, stop going

  ;; effect of toxin on turtle health
  ask turtles[
    if [toxin] of patch-here > 0[ ;; check if toxin is present
      set health (health - (damage-rate * toxin)) ;; update health ;; FOR THE LONG RUN we could reconsider whether this needs another parameter
      if health < death-threshold [ ;; check if health is below threshold
        die ;; if it is, die
      ]]
    if [toxin] of patch-here = 0[ ;; if there's no toxin here
      set health (health + repair-rate) ;; start recovering your health
        if health > 1 [ ;; keep health at or below 1
          set health 1
        ]]]

  ;; toxin degradation
  ask turtles[
    if [toxin] of patch-here > 0 [ ;; if there's toxin on the current patch
      let dr degrade-rate
      ask patch-here [
      set toxin (toxin - dr) ;; degrade some of it ... FOR THE LONG RUN consider whether this should really be linear
        if toxin < 0 [set toxin 0]
  ]]]
  ask patches[ ;; also have to update patch color
    if toxin < toxin-limit [set toxin 0] ;; if toxin goes below threshold, make it zero
    color-patches
  ]

  ;; switching phenotype
  ask turtles[
    if random-float 1 < switch-prob[ ;; choose a random number and compare to switching probability
    set degrade-rate ( ;; determine phenotype (toxin degration capacity) by drawing from the 5-part genotypic distribution
      0.2 * (random-float item 0 sig-vals + (item 0 mu-vals + 0.5 * item 0 sig-vals) +
      random-float item 1 sig-vals + (item 0 mu-vals - 0.5 * item 1 sig-vals) +
      random-float item 2 sig-vals + (item 0 mu-vals - 0.5 * item 2 sig-vals) +
      random-float item 3 sig-vals + (item 0 mu-vals - 0.5 * item 3 sig-vals) +
      random-float item 4 sig-vals + (item 0 mu-vals - 0.5 * item 4 sig-vals)))
      color-turtles
    ]]

  ;; reproduction and mutation
  ask turtles [
    ;; make new cells
    if health > growth-threshold [ ;; only reproduce if your health is good enough
      if (count neighbors with [not any? turtles-here]) > 0 [ ;; only reproduce if there's space
        if random-float 1 < growth-rate [ ;; draw from the probability of reproducing - only if yes, continue
          let growth-space patch-set neighbors with [not any? turtles-here] ;; find patches with no turtles

          ;; store the values that are heritable - for reproduction
          let new-mus [mu-vals] of self
          let new-sigs [sig-vals] of self
          let new-switch-prob [switch-prob] of self
          let new-degrade-rate [degrade-rate] of self
          let new-barcode [barcode] of self

          ;; mutate, if necessary. if mutating, all mus and sigmas may mutate at once. FOR THE LONG RUN we may reconsider that.
          if random-float 1 < mutation-rate [ ; do a random draw to see whether it's time to mutate
            set new-switch-prob (new-switch-prob + ((random-float switch-step) - 0.5  * switch-step)) ; switching probabilty mutates as determined by switch-step
            if new-switch-prob > 1 [ set new-switch-prob 1 ]
            if new-switch-prob < 0 [ set new-switch-prob 0 ]
            ; set new-mus replace-item 0 new-mus ((item 0 new-mus) + ((random-float mu-step) - 0.5 * mu-step))
            foreach [0 1 2 3 4] [ [x] -> ; loop through all mu-values and all sig-values. at each list index between 0 and 4...
              let old-mu item x new-mus ; store the old mu value
              let old-sig item x new-sigs ; store the old sigma value
              let new-mu old-mu + ((random-float mu-step) - 0.5 * mu-step) ; mutate the mu value according to mu-step
              let new-sig old-sig + ((random-float sig-step) - 0.5 * sig-step) ; mutate the sigma value according to sig-step
              if new-mu > (1 - 0.5 * new-sig) [
                set new-mu (1 - 0.5 * new-sig)]
              if new-mu < 0.5 * new-sig[
                set new-mu 0.5 * new-sig]
              set new-mus replace-item x new-mus new-mu ; once you have a good mu-sigma pair, store it
              set new-sigs replace-item x new-sigs new-sig
            ]
          ]

          ;; now reproduce
          ask one-of growth-space [ ;; choose one of the empty patches to sprout
            sprout 1 [ ;; below are all the properties of the newly sprouted cell
              set mu-vals new-mus
              set sig-vals new-sigs
              set switch-prob new-switch-prob
              set degrade-rate new-degrade-rate
              set barcode new-barcode
              set health 1 ;; start with fresh health
              set order-in-stack 1 ;; this will need to be changed if we allow for multiple cells per patch!!
              set growth-rate basal-growth-rate * (1 - tradeoff * degrade-rate ^ alpha)
              ifelse heritable? [ ;; if "heritable" is switched on,
                set degrade-rate new-degrade-rate] ;; inherit degradation rate phenotype from parent
              [set degrade-rate ( ;; otherwise, determine phenotype (toxin degration capacity) by drawing from the 5-part genotypic distribution
                0.2 * (random-float item 0 sig-vals + (item 0 mu-vals - 0.5 * item 0 sig-vals) +
                random-float item 1 sig-vals + (item 0 mu-vals - 0.5 * item 1 sig-vals) +
                random-float item 2 sig-vals + (item 0 mu-vals - 0.5 * item 2 sig-vals) +
                random-float item 3 sig-vals + (item 0 mu-vals - 0.5 * item 3 sig-vals) +
                random-float item 4 sig-vals + (item 0 mu-vals - 0.5 * item 4 sig-vals)))
              ]
              color-turtles
            ]
          ]
        ]
      ]
    ]
  ]


  ;; patches: diffuse toxin + food
  diffuse toxin diff-rate
  ; diffuse food diff-rate ;; un-comment this if we start caring about food
  ask patches[
    if toxin < toxin-limit [set toxin 0] ;; if toxin goes below threshold, make it zero
    color-patches
  ]

  ;; to make toxin get added in periodic pulses
  if pulsing? [
    if ticks-to-pulse <= 0 [ ;; check the ticks-to-pulse counter. if it's at zero...
      pulse ;; add a toxin pulse...
      set ticks-to-pulse pulse-rate ;; and reset the ticks-to-pulse counter
    ]
    set ticks-to-pulse (ticks-to-pulse - 1) ;; if the counter isn't yet at zero, subtract one
  ]

  ;; for turbidostat regime, which dilutes every time the universe fills up
  if turbidostat? [
    if turbidostat-ticks = 20 [ pulse ]
    if count turtles >= count patches [
      dilute
      set turbidostat-ticks 0
    ]
  ]

  set turbidostat-ticks (turbidostat-ticks + 1)
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
0
0
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
11
503
74
536
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

BUTTON
730
47
804
80
poison
poison
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
1096
147
1162
180
dilute
dilute
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
14
408
186
441
dilute-rate
dilute-rate
0
100
5.0
1
1
NIL
HORIZONTAL

SLIDER
729
11
901
44
toxin-conc
toxin-conc
0
1
0.6
.01
1
NIL
HORIZONTAL

BUTTON
813
47
878
80
pulse
pulse
NIL
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
162
104
n
20.0
1
0
Number

PLOT
730
365
930
515
total population
time
cells
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles"

SLIDER
13
108
185
141
basal-growth-rate
basal-growth-rate
0
0.1
0.879
.001
1
NIL
HORIZONTAL

SLIDER
13
145
185
178
damage-rate
damage-rate
0
1
0.7
.01
1
NIL
HORIZONTAL

PLOT
731
213
931
363
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
967
49
1167
199
degradation rate phenotype
degradation rate
NIL
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
180
185
213
tradeoff
tradeoff
0
1
0.71
0.01
1
NIL
HORIZONTAL

SLIDER
12
290
184
323
mu-step
mu-step
0
1
0.51
.01
1
NIL
HORIZONTAL

SLIDER
13
253
185
286
mutation-rate
mutation-rate
0
0.1
0.512
0.001
1
NIL
HORIZONTAL

SLIDER
730
85
902
118
pulse-rate
pulse-rate
0
100
53.0
1
1
NIL
HORIZONTAL

BUTTON
731
121
863
154
start-pulse-rate
start-pulse-rate
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
967
202
1167
352
mu values
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
"mu1" 0.01 0 -955883 true "" "histogram [item 1 mu-vals] of turtles"
"mu0" 0.01 0 -2674135 true "" "histogram [item 0 mu-vals] of turtles"
"mu2" 0.01 0 -13840069 true "" "histogram [item 2 mu-vals] of turtles"
"mu3" 0.01 0 -13345367 true "" "histogram [item 3 mu-vals] of turtles"
"pen-4" 0.01 0 -5825686 true "" "histogram [item 4 mu-vals] of turtles"

PLOT
967
355
1167
505
sigma values
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
"default" 0.01 0 -2674135 true "" "histogram [item 0 sig-vals] of turtles"
"pen-1" 0.01 0 -955883 true "" "histogram [item 1 sig-vals] of turtles"
"pen-2" 0.01 0 -13840069 true "" "histogram [item 2 sig-vals] of turtles"
"pen-3" 0.01 0 -13345367 true "" "histogram [item 3 sig-vals] of turtles"
"pen-4" 0.01 0 -5825686 true "" "histogram [item 4 sig-vals] of turtles"

BUTTON
732
156
863
189
NIL
stop-pulse-rate
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
13
216
185
249
initial-mu
initial-mu
0
1
0.15
.01
1
NIL
HORIZONTAL

SLIDER
12
327
184
360
switch-step
switch-step
0
1
0.2739
0.0001
1
NIL
HORIZONTAL

SWITCH
22
457
154
490
turbidostat?
turbidostat?
1
1
-1000

BUTTON
3
369
66
402
env1
environment-1
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
143
368
206
401
env3
environment-3
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
72
369
135
402
env2
environment-2
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

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
  <experiment name="env3-exp_190629" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="25000"/>
    <metric>count turtles</metric>
    <metric>mean [toxin] of patches</metric>
    <metric>[mu-vals] of turtles</metric>
    <metric>[sig-vals] of turtles</metric>
    <metric>[barcode] of turtles</metric>
    <metric>[degrade-rate] of turtles</metric>
    <enumeratedValueSet variable="toxin-conc">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu-step">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damage-rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="turbidostat?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-mu">
      <value value="0.43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dilute-rate">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pulse-rate">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutation-rate">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="basal-growth-rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tradeoff">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="switch-step">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritable?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sig-step">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-switch-prob">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-sig">
      <value value="0.28"/>
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
