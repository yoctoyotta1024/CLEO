Monoids
=======

Concepts for Templated-Types
--------------------------------
C++ concepts are used abundantly in CLEO. At its simplest, a C++ concept is used to define a set
of constraints on a template. If a certain type satisfies these constraints, it can be called one
possible type of that concept. For example, if a class “Cond” is created to model condensation and
it satisfies the concept of a microphysical process, then “Cond” is a microphysical process.
Likewise another class ‘Colls’ which models collisions is also a microphysical process if it too
satisfies the concept of a microphysical process.

The Use of Concepts to Create Monoids
-------------------------------------
In essence, a monoid is a complete set of members (i.e. elements/types) where each member of the set
has three essential properties:

1. it can be created.

2. it can be combined with another member of the monoid set to produce another member of that
monoid set (i.e. monoid member + monoid member = monoid member).

3. it can be null (a member which does nothing).

As suggested by point (2), there are various different monoids which differ from one another through
their specific definition of these three properties. For example a monoid used to define
"milkshakes" is different to the monoid that defines "melted chocolate" because what it means for (a
member of) the milkshakes to (1)"be created", (2)"be combined" and (3)"be nothing" is not the same
as what it means for (a member of) the melted chocolate.

In CLEO, we make monoids for various different things such as microphysics and observers. We use C++
concepts to ensure only certain types can declare themselves as a member of a particular monoid
(rule(1)); we use some function or structure to define how types which are members of a particular
monoid get combined with one another to produce another member (rule(2)); and we create a null
member of each monoid by making a structure which satisfies the C++ concepts for that monoid but
which does nothing (rule(3)).

We then use our monoids to ensure a templated type does certain things (e.g. microphysics or
observing) and can be combined with other types that also do those things. In this way, we can
create a type, satisfying a given concept, from the combination of several types which also
satisfy that concept and we can use it throughout the code by instantiaing a template
upon compilation. For example, a microphysical process of type ‘CC’ could be created from the
combination of the ‘Cond’ and ‘Colls’ types which are microphysical processes themselves.
Upon compilation, a templated type 'Z' can then be instantiated with 'CC', or 'Cond', or 'Colls'
(or any other member of the microphysics monoid) depending on what we want to model.

A Good Analogy
--------------
The analogy I like to give is mixing paint. Suppose there are a variety of colours;
blue, yellow, red, green, orange etc.. Let’s say that the blue, yellow and red colours
satisfy all the requirements in order for them to be defined as 'wet oil paint' - in other
words, these three colours are wet oil paint.

Meanwhile let's say all the other colours are crayons. The wet oil paints can be mixed in any
combination to create a new wet oil paint - maybe it’s brown, maybe it's violet, maybe it's
something we've never seen before, nevertheless it's certainly wet oil paint. Now of course the red
wet oil paint cannot be mixed with the green crayon to make a new wet oil paint, because clearly
the crayon is not a wet oil paint.

To make the analogy explicit, the requirements used to define wet oil paint and a crayon
are like the C++ concepts used to define a microphysical process or an observer. The
colors are analogous to types satisfying a particular concept, for example ‘Cond’ and
‘Colls’ satisfying the microphysical process. The wet oil paints (and likewise microphysical
processes) not only satisfy a concept but also are monoids. They therefore have a specified rule
(function) which allows them to be combined to create new wet oil paints, analogous to
creating the microphysical process‘CC’ from ‘Cond’ and ‘Colls’. (The null property of a monoid
would be like transparent wet oil paint - it can be used like all wet oil paints, but it
doesn’t do anything.)
