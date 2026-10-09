{{ name | escape | underline}}

.. automodule:: {{ fullname }}
   :members:

{# autodoc cannot see argparse options, so render the parser itself #}
{% if fullname == 'gfnff.cli' %}
Options
-------

.. argparse::
   :module: gfnff.cli
   :func: _build_parser
   :prog: gfnff
   :nodescription:
{% endif %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: module_template.rst
   :recursive:
{# pkgutil lists the bundled libgfnff.so as a submodule; it is not one #}
{% for item in modules if item.split('.')[-1] != 'libgfnff' %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
