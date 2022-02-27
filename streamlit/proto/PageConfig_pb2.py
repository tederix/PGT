# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: streamlit/proto/PageConfig.proto

from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor.FileDescriptor(
  name='streamlit/proto/PageConfig.proto',
  package='',
  syntax='proto3',
  serialized_options=None,
  create_key=_descriptor._internal_create_key,
  serialized_pb=b'\n streamlit/proto/PageConfig.proto\"\x97\x03\n\nPageConfig\x12\r\n\x05title\x18\x01 \x01(\t\x12\x0f\n\x07\x66\x61vicon\x18\x02 \x01(\t\x12\"\n\x06layout\x18\x03 \x01(\x0e\x32\x12.PageConfig.Layout\x12\x37\n\x15initial_sidebar_state\x18\x04 \x01(\x0e\x32\x18.PageConfig.SidebarState\x12)\n\nmenu_items\x18\x05 \x01(\x0b\x32\x15.PageConfig.MenuItems\x1a\x87\x01\n\tMenuItems\x12\x14\n\x0cget_help_url\x18\x01 \x01(\t\x12\x15\n\rhide_get_help\x18\x02 \x01(\x08\x12\x18\n\x10report_a_bug_url\x18\x03 \x01(\t\x12\x19\n\x11hide_report_a_bug\x18\x04 \x01(\x08\x12\x18\n\x10\x61\x62out_section_md\x18\x05 \x01(\t\" \n\x06Layout\x12\x0c\n\x08\x43\x45NTERED\x10\x00\x12\x08\n\x04WIDE\x10\x01\"5\n\x0cSidebarState\x12\x08\n\x04\x41UTO\x10\x00\x12\x0c\n\x08\x45XPANDED\x10\x01\x12\r\n\tCOLLAPSED\x10\x02\x62\x06proto3'
)



_PAGECONFIG_LAYOUT = _descriptor.EnumDescriptor(
  name='Layout',
  full_name='PageConfig.Layout',
  filename=None,
  file=DESCRIPTOR,
  create_key=_descriptor._internal_create_key,
  values=[
    _descriptor.EnumValueDescriptor(
      name='CENTERED', index=0, number=0,
      serialized_options=None,
      type=None,
      create_key=_descriptor._internal_create_key),
    _descriptor.EnumValueDescriptor(
      name='WIDE', index=1, number=1,
      serialized_options=None,
      type=None,
      create_key=_descriptor._internal_create_key),
  ],
  containing_type=None,
  serialized_options=None,
  serialized_start=357,
  serialized_end=389,
)
_sym_db.RegisterEnumDescriptor(_PAGECONFIG_LAYOUT)

_PAGECONFIG_SIDEBARSTATE = _descriptor.EnumDescriptor(
  name='SidebarState',
  full_name='PageConfig.SidebarState',
  filename=None,
  file=DESCRIPTOR,
  create_key=_descriptor._internal_create_key,
  values=[
    _descriptor.EnumValueDescriptor(
      name='AUTO', index=0, number=0,
      serialized_options=None,
      type=None,
      create_key=_descriptor._internal_create_key),
    _descriptor.EnumValueDescriptor(
      name='EXPANDED', index=1, number=1,
      serialized_options=None,
      type=None,
      create_key=_descriptor._internal_create_key),
    _descriptor.EnumValueDescriptor(
      name='COLLAPSED', index=2, number=2,
      serialized_options=None,
      type=None,
      create_key=_descriptor._internal_create_key),
  ],
  containing_type=None,
  serialized_options=None,
  serialized_start=391,
  serialized_end=444,
)
_sym_db.RegisterEnumDescriptor(_PAGECONFIG_SIDEBARSTATE)


_PAGECONFIG_MENUITEMS = _descriptor.Descriptor(
  name='MenuItems',
  full_name='PageConfig.MenuItems',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='get_help_url', full_name='PageConfig.MenuItems.get_help_url', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=b"".decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='hide_get_help', full_name='PageConfig.MenuItems.hide_get_help', index=1,
      number=2, type=8, cpp_type=7, label=1,
      has_default_value=False, default_value=False,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='report_a_bug_url', full_name='PageConfig.MenuItems.report_a_bug_url', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=b"".decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='hide_report_a_bug', full_name='PageConfig.MenuItems.hide_report_a_bug', index=3,
      number=4, type=8, cpp_type=7, label=1,
      has_default_value=False, default_value=False,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='about_section_md', full_name='PageConfig.MenuItems.about_section_md', index=4,
      number=5, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=b"".decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=220,
  serialized_end=355,
)

_PAGECONFIG = _descriptor.Descriptor(
  name='PageConfig',
  full_name='PageConfig',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='title', full_name='PageConfig.title', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=b"".decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='favicon', full_name='PageConfig.favicon', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=b"".decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='layout', full_name='PageConfig.layout', index=2,
      number=3, type=14, cpp_type=8, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='initial_sidebar_state', full_name='PageConfig.initial_sidebar_state', index=3,
      number=4, type=14, cpp_type=8, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='menu_items', full_name='PageConfig.menu_items', index=4,
      number=5, type=11, cpp_type=10, label=1,
      has_default_value=False, default_value=None,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[_PAGECONFIG_MENUITEMS, ],
  enum_types=[
    _PAGECONFIG_LAYOUT,
    _PAGECONFIG_SIDEBARSTATE,
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=37,
  serialized_end=444,
)

_PAGECONFIG_MENUITEMS.containing_type = _PAGECONFIG
_PAGECONFIG.fields_by_name['layout'].enum_type = _PAGECONFIG_LAYOUT
_PAGECONFIG.fields_by_name['initial_sidebar_state'].enum_type = _PAGECONFIG_SIDEBARSTATE
_PAGECONFIG.fields_by_name['menu_items'].message_type = _PAGECONFIG_MENUITEMS
_PAGECONFIG_LAYOUT.containing_type = _PAGECONFIG
_PAGECONFIG_SIDEBARSTATE.containing_type = _PAGECONFIG
DESCRIPTOR.message_types_by_name['PageConfig'] = _PAGECONFIG
_sym_db.RegisterFileDescriptor(DESCRIPTOR)

PageConfig = _reflection.GeneratedProtocolMessageType('PageConfig', (_message.Message,), {

  'MenuItems' : _reflection.GeneratedProtocolMessageType('MenuItems', (_message.Message,), {
    'DESCRIPTOR' : _PAGECONFIG_MENUITEMS,
    '__module__' : 'streamlit.proto.PageConfig_pb2'
    # @@protoc_insertion_point(class_scope:PageConfig.MenuItems)
    })
  ,
  'DESCRIPTOR' : _PAGECONFIG,
  '__module__' : 'streamlit.proto.PageConfig_pb2'
  # @@protoc_insertion_point(class_scope:PageConfig)
  })
_sym_db.RegisterMessage(PageConfig)
_sym_db.RegisterMessage(PageConfig.MenuItems)


# @@protoc_insertion_point(module_scope)
