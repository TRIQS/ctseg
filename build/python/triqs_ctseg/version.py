################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2014-2018, P. Seth, I. Krivenko, M. Ferrero, O. Parcollet
# Copyright (C) 2018-2019, Simons Foundation
#   author: N. Wentzell
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

version = "3.1.0"
triqs_hash = "380bf41014b53651f4cd339f4242c98bc9048b0e"
triqs_ctseg_hash = "6983f08f8c4f0504fceacf77adc2c4f81fce33a4"

def show_version():
  print("\nYou are using triqs_ctseg version %s\n"%version)

def show_git_hash():
  print("\nYou are using triqs_ctseg git hash %s based on triqs git hash %s\n"%("6983f08f8c4f0504fceacf77adc2c4f81fce33a4", triqs_hash))
