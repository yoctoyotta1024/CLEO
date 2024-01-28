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
Within CLEO, a monoid can be thought of as a type that has three essential properties:

1. it can be created.

2. it can be combined with a monoid of its kind to produce a new
monoid of that kind (i.e. a monoid + a monoid = a monoid).

3. it can be null (a structure which does nothing).

As suggested by point (2), there are various kinds of monoids. The different kinds of monoids
differ over the specifics of these three properties. For example the kind of monoid which enacts
microphysics is different to the kind of monoid that does observation because what it means for
a microphysics monoid to (1)"be created", (2)"be combined" and (3)"be nothing" is not the same
as what it means for an observation monoid.

For each kind of monoid, we use C++ concepts when we want to make certain a templated type is
of that kind. In these situations, the rule for combining that kind of monoid is implemented by
some function and the concepts are used to ensure that only types which obey that rule are allowed
to declare themselves as that kind of monoid. By doing this, we can then create a monoid,
satisfying a given concept, from the combination of several types which also
satisfy that concept. For example, a microphysical process of type ‘CC’ could be created from the
combination of the ‘Cond’ and ‘Colls’ types which are microphysical processes themselves.

A Good Analogy
--------------
The analogy I like to give is mixing paint. Suppose there are a variety of colours;
blue, yellow, red, green, orange etc.. Let’s say that the blue, yellow and red colours
satisfy all the requirements in order for them to be defined as 'wet oil paint' - in other
words, these three colours are wet oil paint. Meanwhile let's say all the other colours
are crayons. The wet oil paints can be mixed in any combination to create a new wet oil paint -
maybe it’s brown, maybe it's violet, maybe it's something we've never seen before, nevertheless
it's certainly wet oil paint. Now of course the red wet oil paint cannot be mixed with the
green crayon to make a new wet oil paint, because clearly the crayon is not a wet oil paint. To
make the analogy explicit, the requirements used to define wet oil paint and a crayon
are like the C++ concepts used to define a microphysical process or an observer. The
colors are analogous to types satisfying a particular concept, for example ‘Cond’ and
‘Colls’ satisfying the microphysical process. The wet oil paints (and likewise microphysical
processes) not only satisfy a concept but also are monoids. They therefore have a specified rule
(function) which allows them to be combined to create new wet oil paints, analogous to
creating the microphysical process‘CC’ from ‘Cond’ and ‘Colls’. (The null property of a monoid
would be like transparent wet oil paint - it can be used like all wet oil paints, but it
doesn’t do anything.)
